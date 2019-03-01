function [J,I,V] = comp_u_tri(U)
%Find the compressed sparse column format of a upper triangular matrix
%   Input  - U  upper triangular matrix
%   Output - J  row indices
%            I  column pointers
%            V  nonzero entries
[J,~,V] = find(U); I = ones(size(U,1)+1,1);
for k = 2:length(I); I(k) = I(k-1) + nnz(U(:,k-1)); end
end