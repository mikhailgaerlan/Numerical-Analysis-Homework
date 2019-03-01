function [J,I,V] = comp_l_tri(L)
%Find the compressed sparse column format of a lower triangular matrix
%   Input  - L  lower triangular matrix
%   Output - J  row indices
%            I  column pointers
%            V  nonzero entries
SL = tril(L,-1);
[J,~,V] = find(SL); I = ones(size(SL,1),1);
for k = 2:length(I); I(k) = I(k-1) + nnz(SL(:,k-1)); end
end