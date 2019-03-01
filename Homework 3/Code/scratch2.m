m = 99;
h = 1/(m+1);
[X,Y]=meshgrid(h:h:1-h);
x = X(1,:);y = Y(:,1);
F = f(X,Y);
b0 = v(x,0);b1 = v(x,1);
c0 = v(0,y);c1 = v(1,y);
Vapprox = solvePoisson(m,F,b0,b1,c0,c1);
Vactual = v(X,Y);
max(abs(Vapprox(:)-Vactual(:)))
colormap([1 0 0; 0 0 1;0 1 0]);
figure
s = surf(X,Y,Vapprox);hold on
t = surf(X,Y,Vactual);hold on
u = surf(X,Y,F);

function true = v(x,y)
true = x.^3.*y.^2;
end

function fun = f(x,y)
fun = -6.*x.^2.*y.^2-2.*x.*y.^3;
end