%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Lagrange Interpolation Function
% 
% Author: Chad McKell
% Date: 14 Dec 2016
% Place: University of Edinburgh
%
% Description: This function produces a matrix of QxN Langrange polynomials.
% Each of the N constituent polynomials of the global interpolation
% polynomial is sampled over a range of Q discrete values. The matrix is 
% suitable for performing Nth-order Lagrange interpolation.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function L = interpLG(Q,N)
% Q: number of discrete values over which the polynomials are sampled.
% N: number of Lagrange polynomials used for interpolation. Each polynomial
% corresponds to a separate column of the matrix L. The highest order 
% polynomial is N-1 order.

% Vector of polynomial estimation points
a = (-(N-1):2:N-1)/2;

% Range limit over which polynomials are evaluated
range = 0.5;

% Vector of bins spanning interpolation range. Each bin is separated by 1/Q.
alpha = -range:1/Q:range;

% The last of value of alpha is removed in order to define the correct
% interpolation range of [-range:range-(1/Q]
alpha = alpha(1:end-1);

% Initialized Lagrange polynomial matrix with Q rows and N columns
L = ones(Q,N);

% For each bin q along the interpolation range, calculate the N Langrange
% polynomials. Each polynomial is given by the product of the expression
% (alpha(q)-a(k))/(a(n)-a(k)) with itself for every n is not equal to k,
% where n and k index the vector 'a' of estimation points. 
for q = 1:Q
    for n = 1:N
        for k = 1:N
            if n ~= k
                L(q,n) = L(q,n)*(alpha(q)-a(k))/(a(n)-a(k));
            end
        end
    end
end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Comments
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Reference 3 was a very helpful guide to constructing the matrix in
% this MATLAB assignment. However, since the polynomials are sampled at Q 
% discrete values in the [-range:range-(1/Q], the matrix had to interate 
% from 1:Q.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% References
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% 1. Newton, M. MAFTDSP Course Notes, Lecture 10: Systems, the DTFT, Simple 
% Digital Filters & Interpolation, University of Edinburgh (2016).
% 2. Website: http://www.math.usm.edu/lambers/mat772/fall10/lecture5
% 3. Website: http://uk.mathworks.com/matlabcentral/fileexchange/899-...
% ...lagrange-polynomial-interpolation


