function y = toeplitz(t,x)
%Multiply a column vector by a Toeplitz matrix using fft.
%   input  - t    row vector that generates a Toeplitz matrix, T
%          - x    column vector
%   output - T*x  matrix multiplication
n = length(x);
c = zeros(2*n-1,1);
c(1:n) = t(n:end);
c(n+1:end) = t(1:n-1);

d = zeros(2*n-1,1);
d(1:n) = x;
y = circulant(c',d);
y = y(1:n);
end

