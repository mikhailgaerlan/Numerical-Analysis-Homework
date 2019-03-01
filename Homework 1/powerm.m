function [lambda,j,relresid,x] = powerm(n,Q,Jv,x0,eps,kmax)
%Find the dominant eigenvalue of a matrix A'
%   input  - A         n x n matrix
%          - x0        starting vector
%          - eps       tolerance
%          - kmax      maximum iterations
%   output - lambda    dominant eigenvalue of A'
%          - j         iteration reached
%          - relresid  relative residual
x = x0;
[~,I]=max(abs(x));lambda = x0(I);
for k = 1:kmax
    j = k; oldlambda = lambda;
    x = graphx(n,Q,Jv,x);
    [~,I]=max(abs(x));
    lambda=x(I); x = x/lambda;
    relresid=abs(lambda-oldlambda)/abs(oldlambda);
    if (relresid <= eps)
        break;
    end
end
end