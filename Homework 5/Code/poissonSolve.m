function [v,iter,relres] = poissonSolve(m,b,gamma,tol)
%Solve two-dimensional Poisson's equation with extra x-derivative
%   input  - m            matrix size
%            b            right-hand side vector
%            gamma        parameter
%            tol          error tolerance
%            v0           initial vector
%   output - v            solution vector
%          - iter         number of iterations
%          - relres       relative residual norm of final GMRES iterate
[v,~,relres,iters,~] = gmres(@(x)multAp(x,gamma,m),solveM1(b),[],tol,m^2);
iter = iters(2);
    function v = solveM1(x)
        H = reshape(x,[m,m]);
        V = solvePoisson(m,H);
        v = reshape(V,[m*m,1]);
    end
end
