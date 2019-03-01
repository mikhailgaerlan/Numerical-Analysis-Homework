function x = toeplitzpcg(t,b)
%Solve a system Tx=b where T is symmetric positive-definite Toeplitz.
%   input  - t      row vector for symetric Toeplitz matrix, T
%          - b      column vector
%          - tol    tolerance
%          - maxit  maximum iterations
%   output - x      vector
    x = pcg(@toeplitzm,b,1e-9,10000);
    
    function y = toeplitzm(x)
        y = toeplitz([t(end:-1:2),t],x);
    end
end