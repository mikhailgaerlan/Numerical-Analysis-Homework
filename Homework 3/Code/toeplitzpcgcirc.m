function x = toeplitzpcgcirc(t,b)
%Solve a system Tx=b where T is symmetric positive-definite Toeplitz.
%Uses a right preconditioner C^-1.
%   input  - t      row vector for symetric Toeplitz matrix, T
%          - b      column vector
%          - tol    tolerance
%   output - x      solution
n = length(b); c = ((n:-1:1).*t+(0:(n-1)).*t([1,end:-1:2]))/n;
x = pcg(@toeplitzm,b,1e-9,10000,@circulantinv);

    function y = circulantinv(xx)
        y = conj(fft(conj(((conj(fft(c'))).^(-1)).*fft(xx))))/n;
    end

    function y = toeplitzm(xx)
        y = toeplitz([t(end:-1:2),t],xx);
    end
end