format long e; clearvars; clc; m = 7; nmax = 100;

params = [1/2,0,0,1,1/2;...
    1/2,0,2,7/2,2;...
    1/4,1,0,1,1/2;...
    1/4,1,2,7/2,2;...
    1/2,2,0,1,1/2;
    1/2,2,2,7/2,2];
n = size(params,1);
tols = []; ks = [];

for i = 1:n
    a=params(i,1);b=params(i,2);
    alpha=params(i,3);beta=params(i,4);gamma=params(i,5);
    [X,Y,Vactual,Vapprox,tol,k] = mult_Schwarz(m,nmax,a,b,@(x,y) f(x,y,alpha,beta,gamma),@(x,y) v(x,y,alpha,beta,gamma));
    tols = [tols;tol]; ks = [ks;k];
    figure(1);scatter3(X,Y,Vactual,'square','CData',Vactual);view(2);
    figure(2);scatter3(X,Y,Vapprox,'square','CData',Vapprox);view(2);
end
tols
ks

function true = v(x,y,a,b,g)
true = y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end

function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end