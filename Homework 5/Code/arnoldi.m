function [Hkt,Vkk] = arnoldi(A,r,kmax,tol)
%Use the Arnoldi process on a matrix A
%   input  - A     matrix
%            r     starting vector
%            kmax  max iterations
%            tol   norm(q) tolerance
%   output - Hkt   upper-Hessenberg matrix
%            Vkk   matrix of Arnoldi vectors
beta = norm(r); Vkk = r/beta;
for k = 1:kmax
    q = A(Vkk(:,k));
    for j = 1:k; Hkt(j,k) = Vkk(:,j)'*q; q = q-Hkt(j,k)*Vkk(:,j); end
    Hkt(k+1,k) = norm(q);
    if norm(q) <= tol; break; end
    Vkk(:,k+1) = q/Hkt(k+1,k);
end
end

