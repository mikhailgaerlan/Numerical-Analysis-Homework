function x = solve_lower(c,J,I,V)
%Solve Lx = c where L is a lower-triangular matrix
%   Input  - J  row indices
%            I  pointers
%            V  nonzero entries
%            c  right-hand side
%   Output - x  solution
x = c;
for k = 1:length(x)
    x(k) = x(k)/V(I(k));
    indices = (I(k)+1):(I(k+1)-1); rows = J(indices);
    x(rows) = x(rows) - x(k)*V(indices);
end
end