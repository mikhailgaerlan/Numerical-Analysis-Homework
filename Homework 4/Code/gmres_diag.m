function [x,m,relres] = gmres_diag(A,b,k0,tol,x0)
%Solve Ax=b using GMRES with a right diagonal preconditioner
%   Input  - A       matrix
%          - b       right-hand side
%          - k0      restart parameter
%          - tol     tolerance
%          - x0      initial solution
%   Output - x       approximate solution
%          - m       number of matrix-vector products
%          - relres  history of relative residual norms
[x,~,~,iter,relres] = gmres(A,b,k0,tol,length(b),@idens,@diags,x0);
m = (iter(1)-1)*k0+iter(2)+iter(1);relres = relres/norm(b-A*x0);
    function x = idens(b); x = b; end
    function x = diags(b); x = b./diag(A); end
end