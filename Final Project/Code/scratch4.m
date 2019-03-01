clearvars; clf;
a=0;b=1;c=0;d=5;
alpha = 2;beta = 7/2;gamma = 2;
%m1 = 15; m2 = 3;
%h1 = (b-a)/(m1+1); h2 = (d-c)/(m2+1);
m = 10; h1 = 1/2^m; h2 = 1/2^m;

[X,Y]=meshgrid((a+h1):h1:(b-h1),(c+h2):h2:(d-h2));
Vapprox = poissonRectSolve(a,b,c,d,h1,h2,@(x,y) f(x,y,alpha,beta,gamma),@(x,y) v(x,y,alpha,beta,gamma));

Vactual = v(X,Y,alpha,beta,gamma);
figure(1); s = surf(X,Y,Vactual); set(s,'LineStyle','none');
figure(2); s = surf(X,Y,Vapprox); set(s,'LineStyle','none');
max(abs(Vapprox(:)-Vactual(:)))/max(abs(Vactual(:)))

function true = v(x,y,a,b,g)
true = y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end

function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end