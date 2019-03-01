function [J,I] = get_lower(A)
%Find the sparse column format of a the lower-triangular part of a matrix
%   Input  - A  sparse matrix
%   Output - J  row indices
%            I  column pointers
%            V  nonzero entries
[J,K] = find(tril(A)); I = find(J-K==0);I = [I;nnz(tril(A))+1];
end