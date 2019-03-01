function V = iCholesky(A,J,I)
%Find the elements of an incomplete Cholesky matrix
%   Input  - A  sparse matrix
%            J  row indices
%            I  column pointers
%   Output - V  entries of incomplete Cholesky factorization
n = size(A,1); [~,~,V] = find(tril(A));
for k = 1:n
    if V(I(k)) <= 0; break; end
    V(I(k))=sqrt(V(I(k)));
    indJ = (I(k)+1):(I(k+1)-1);
    V(indJ) = V(indJ)/V(I(k));
    for j = indJ
        indI = I(J(j)):(I(J(j)+1)-1);
        [~,rowsJ,rowsI] = intersect(J(indJ),J(indI));
        rowsJ = rowsJ + I(k); rowsI = rowsI + I(J(j)) - 1;
        vj = V(j); V(rowsI) = V(rowsI) - vj*V(rowsJ);
    end
end
end