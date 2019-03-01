function c = solve_lowert(b,J,I,V)
%Solve L^Tc=b where L is lower-triangular
%   Input  - J  row indices
%            I  pointers
%            V  nonzero entries
%            b  right-hand side
%   Output - c  solution
c = b;c(end) = c(end)/V(end);
for k = (length(c)-1):-1:1
    columns = (I(k)+1):(I(k+1)-1); rows = J(columns);
    c(k) = (c(k) - sum(V(columns).*c(rows)))/V(I(k));
end
end