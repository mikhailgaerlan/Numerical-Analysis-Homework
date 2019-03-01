function x = solve_u(J,I,V,c)
%Solve Ux = c where U is a nonsingular upper-triangular
%   Input  - J  row indices
%            I  pointers
%            V  nonzero entries
%            c  right-hand side
%   Output - x  solution
x = c;
for k = (length(I)-1):-1:2
    indices = I(k):(I(k+1)-1);rows = J(indices);
    x(k) = x(k)/V(indices(end));
    x(rows(1:end-1)) = x(rows(1:end-1)) - V(indices(1:end-1))*x(k);
end; x(1) = x(1)/V(1);
end