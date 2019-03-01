function y = circulant(c,x)
%Multiply a column vector by a circulant matrix using fft.
%   input  - c    row vector that generates a circulant matrix, C
%          - x    column vector
%   output - C*x  matrix multiplication

n = length(x);
y = conj(fft(conj(conj(fft(c')).*fft(x))))/n;
end