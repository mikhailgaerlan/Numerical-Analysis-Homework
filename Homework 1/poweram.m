function [lambda,j,relresid,x] = poweram(n,Q,Jv,a,x0,eps,kmax)
%Find the dominant eigenvalue of a matrix A_a'
%   input  - A         n x n matrix
%          - x0        starting vector
%          - eps       tolerance
%          - kmax      maximum iterations
%   output - lambda    dominant eigenvalue of A_a'
%          - j         iteration reached
%          - relresid  relative residual
x = x0;
[~,I]=max(abs(x));lambda = x0(I);
for k = 1:kmax
    j = k; oldlambda = lambda;
    x = graphax(n,Q,Jv,a,x);
    [~,I]=max(abs(x));
    lambda=x(I); x = x/lambda;
    relresid=abs(lambda-oldlambda)/abs(oldlambda);
    if (relresid <= eps)
        break;
    end
end
end