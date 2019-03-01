function Vapprox = poissonRectSolve(a,b,c,d,h1,h2,f,v)
%Solve Poisson equation -u_xx-u_yy=f(x,y)
%   input  - f            right-hand side function
%            v            boundary conditions
%            [a,b]        x interval
%            [c,d]        y interval
%            [h1,h2]      x and y intervals
%   output - Vapprox      solution matrix

[X,Y]=meshgrid((a+h1):h1:(b-h1),(c+h2):h2:(d-h2));
x = X(1,:);y = Y(:,1); F = f(X,Y);
b0 = v(x,c);c0 = v(a,y);b1 = v(x,d);c1 = v(b,y);
F = h1^2*h2^2*F;
F(1,:) = F(1,:) + h1^2*b0; F(end,:) = F(end,:) + h1^2*b1;
F(:,1) = F(:,1) + h2^2*c0; F(:,end) = F(:,end) + h2^2*c1;

Vapprox = solvePoissonRect(F,a,b,c,d);

end