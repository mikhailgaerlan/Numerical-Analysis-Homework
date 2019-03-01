function c = solve_l(J,I,V,b)
%Solve Lc = b where L is unit lower-triangular
%   Input  - J  row indices
%            I  pointers
%            V  nonzero entries
%            b  right-hand side
%   Output - c  solution
c = b;
for k = 1:(length(I)-1)
    indices = I(k):(I(k+1)-1); rows = J(indices);
    c(rows) = c(rows) - V(indices)*c(k);
end
end