function [x,m,relres] = gmres_ssor(A,b,k0,tol,D,x0)
%Solve A'x'=b' using GMRES with an SSOR-type preconditioner
%   Input  - A       matrix
%          - b       right-hand side
%          - tol     tolerance
%          - k0      restart parameter
%          - D       diagonal matrix
%          - x0      initial solution
%   Output - x       approximate solution
%          - m       number of matrix-vector products
%          - relres  history of relative residual norms
D0 = diag(diag(A)); D1 = D0 - 2*D;
F = -tril(A,-1); G = -triu(A,1);M1 = (D-F)*(D^(-1)); M2 = D-G; 
[x,~,~,iter,relres] = gmres(@mult_Ap,M1\b,k0,tol,length(b),[],[],M2*x0);
m = (iter(1)-1)*k0+iter(2)+iter(1);relres = relres/norm(b-A*x0);
    function b = mult_Ap(x)
        u = (D-G)\x;
        l = (D-F)\(x+diag(D1).*u);
        b = diag(D).*(u+l);
    end
end