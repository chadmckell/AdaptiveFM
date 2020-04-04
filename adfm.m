%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Adaptive FM Synthesis
% 
% Author: Chad McKell
% Date: 10 Dec 2019
% Place: University of California San Diego
%
% Description: This script modifies the harmonic content of a signal using
% adaptive FM synthesis.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function adfm(input)
% input: type of audio input ('flute' or 'clarinet')
tic;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Read in the input .WAV file 'x' and sample rate 'Fs'
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
switch input
    
    case 'flute'
        [x,Fs] = audioread('Flute_nonvib_ff_C4.wav'); 
        f0 = 261; % fundamental frequency (Hz)

    case 'clarinet'
        [x,Fs] = audioread('BbClar_ff_D3.wav'); 
        f0 = 146.5; % fundamental frequency (Hz)
end

% If user inserts stereo audio, only read the left channel 
if size(x,2) > 1
    x = x(:,1); 
    warning('File x has two channels. Only the left channel was read.');
end

% Find input sample length
L = length(x);  

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Define global parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fm = f0;                % modulation frequency (Hz)
fc = f0;                % carrier frequency (Hz)
g = 0.5;                % effect strength, set between -1 and 1
Ts = 1/Fs;              % sampling period
n = 1:L;                % sample bin vector              

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Define variable delay
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Imax = 1;                                       % max index of modulation
Ienv = Imax * adsr(0.05, 0.05, 0.1, 0.8, L);    % index of modulation
beta = Ienv * 0.5*Fs/pi/fc;                     % multiplicative factor
D = beta + beta.*cos(2*pi*fm*n*Ts);             % variable delay

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Generate Lagrange interpolation matrix and evaluation parameters. (Note
% that the delay-line and interpolation code below comes directly from  
% previous work I completed as a student at the University of Edinburgh).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Set Lagrange interpolation parameters   
Q = 100; % number of discrete values over which polynomials are sampled
N = 4; % number of polynomials used for interpolation

LG = interpLG(Q,N);

% Define range limit over which polynomials are evaluated 
range = 0.5; 

% Define vector of evaluation points spanning interpolation range
alpha = -range:1/Q:range; 

% Remove the last value of alpha in order to define the correct
% interpolation range of [-range:range-(1/Q)]
alpha = alpha(1:end-1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Produce modified output signal
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize output signal
y = zeros(L,1); 

% Find max value of oscillator D. Convert it to an integer.
maxD = round(max(D), 0);

% Let y equal x from n=1 to n=maxD
y(1:maxD) = x(1:maxD);

% Compute the contribution to yFL from the LFO by interpolating the value
% of x(n-D(n)), multiplying it by g, then adding it to x(n).
for n = maxD+1:L-1
    
    % Compute the closest integer value BELOW the total index 'n-D(n)'
    Dint = floor(n-D(n));
    
    % Compute the closest integer value ABOVE the total index 'n-D(n)'
    Dint1 = ceil(n-D(n));

    % Construct vector of polynomial estimation points centered around the
    % midway point between Dint and Dint1
    a = Dint+0.5 + (-(N-1):2:N-1)/2;
    
    % Construct vector of polynomial coefficients evaluated at the 
    % polynomial estimation points defined above
    x_alpha = x(a);

    % Locate the point along n to be interpolated (as defined in the alpha
    % coordinate system centered around the midway point between the values
    % Dint and Dint1). 
    alphaD = (n-D(n))-Dint-0.5;

    % Find the index I in the vector 'alpha' that is nearest to alphaD. 
    [~,I] = min(abs(alpha(1:end) - alphaD));

    % Find the estimate of x(n-D(n)) by taking the dot product of the Ith 
    % row in the Langrange polynomial matrix 'LG' with the vector 'x_alpha'
    % of polynomial coefficients 
    xD = LG(I,:)*x_alpha;
     
    % Multiply the estimate of x(n-D(n)) by the gain parameter g.
    y(n) = g*xD;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Listen to output signal 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soundsc(y,Fs); 


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Save simulation as .wav file
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define filename. Include parameter values in filename.
filename = 'adfm.wav';

% Write .wav file to MATLAB pwd at 16 bits
audiowrite(filename, y*2, Fs, 'BitsPerSample', 16);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plot original and modified signals
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
plotspec(y, Fs, 'db');
hold on;
plotspec(x, Fs, 'db');
lgd = legend('Modified', 'Original');
lgd.FontSize = 12;
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Check code efficiency
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
toc % print elapsed time

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% References
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% 1. Newton, M. MAFTDSP Course Notes, Lecture 10: Systems, the DTFT, Simple 
% Digital Filters & Interpolation, University of Edinburgh (2016).
% 2. Newton, M. MAFTDSP Course Notes, Tutorial 8: Windows, the DFT, and 
% Analysis-Synthesis, University of Edinburgh (2016).

end


function env = adsr(a,d,r,sus,L)

% a: length of attack (expressed as ratio of complete signal)
% d: length of decay (expressed as ratio of complete signal)
% r: length of release (expressed as ratio of complete signal)
% sus: sustain level (value between 0 and 1)
% N: number of samples in signal

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute duration of each section of the envelope in samples
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

a_dur = round(L*a);                           % attack duration (samples)
d_dur = round(L*d);                           % decay duration (samples)
r_dur = round(L*r);                           % release duration (samples)
s_dur = L - a_dur - d_dur - r_dur;            % sustain duration (samples)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute envelope
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

att = 0:1/a_dur:1;                            % attack function
dec = 1:(sus-1)/d_dur:sus;                    % decay function
rel = sus:-sus/r_dur:0;                       % release function
sus = sus*ones(1, s_dur-3);                   % sustain function (delete 3 extra samples)
env = [ att, dec, sus, rel ];                 % complete envelope

end


function plotspec(x, fs, option)
%
% PLOTSPEC Plot the magnitude and phase of the fft.
%
% X is the input signal
% FS is the sampling rate
% 'OPTION' determines on which scale (linear or dB) the magnitude should be
% plotted.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the DFT of the signal x
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = length(x);
y = fft(x, N);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the magnitude and phase of the spectrum y
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y_mag = abs(y);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Normalize yMag
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y_mag = y_mag/max(y_mag);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the frequency axis (note that this only works if the sampling 
% rate you passed in to plotspec matches the sampling rate you used to  
% define the input signal x).
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f = [0:N-1]*fs/N;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot the magnitude and phase of the spectrum y
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(1);
if strcmp(option, 'linear')
    
    plot(f, y_mag, 'LineWidth', 1.5);
    xlim([0,15000]);
    grid on;
    title('Magnitude Response');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (normalized)');
    
elseif strcmp(option,'db')
    
    plot(f, 10*log10(y_mag), 'LineWidth', 1.5); 
    xlim([0,15000]); 
    grid on;
    title('Magnitude Response');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
else
    error('You must choose linear or db for the option');
end

end
 
