format long e; clearvars; clc; m = 4; nmax = 100; tol = 1e-10;
a = 1/2; b = 2; alpha = 2; beta = 7/2; gamma = 2;

[X,Y,Vactual,Vapprox,~,~] = mult_Schwarz(m,nmax,a,b,@(x,y) f(x,y,alpha,beta,gamma),@(x,y) v(x,y,alpha,beta,gamma),tol);

[x,y] = meshgrid(X,Y);
figure(1); surf(X,Y,Vactual)

function true = v(x,y,a,b,g)
true = y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end
tmatla
function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end