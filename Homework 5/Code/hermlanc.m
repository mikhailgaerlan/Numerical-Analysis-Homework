function Tk = hermlanc(A,r,kmax,tol)
%Apply the Hermitian Lanczos process for a Hermitian matrix
%   input  - A     Hermitian matrix
%          - r     starting vector
%          - kmax  max iterations
%          - tol   norm(q) tolerance
%   output - Tk    tridiagonal matrix
beta = norm(r); vk = r/beta;
for k = 1:kmax
    q = A*vk; if k > 1; q = q - beta*vk1; end
    Tk(k,k) = vk'*q; q = q - Tk(k,k)*vk; beta = norm(q);
    if beta <= tol; break; end
    if k ~= kmax
        Tk(k+1,k) = beta; Tk(k,k+1) = beta;
        vk1 = vk; vk = q/beta;
    end
end
end

