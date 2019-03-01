clearvars; clf;
c=0;d=6;a=0;b=8;
alpha = 1;beta = 2;gamma = 1/2;
m1 = 500; m2 = 800;
h1 = (b-a)/(m1+1); h2 = (d-c)/(m2+1);

[X,Y]=meshgrid((a+h1):h1:(b-h1),(c+h2):h2:(d-h2));
x = X(1,:);y = Y(:,1);
F = f(X,Y,alpha,beta,gamma);
b0 = v(x,c,alpha,beta,gamma);c0 = v(a,y,alpha,beta,gamma);
b1 = v(x,d,alpha,beta,gamma);c1 = v(b,y,alpha,beta,gamma);
F = h1^2*h2^2*F;
F(1,:) = F(1,:) + h1^2*b0; F(end,:) = F(end,:) + h1^2*b1;
F(:,1) = F(:,1) + h2^2*c0; F(:,end) = F(:,end) + h2^2*c1;

Vactual = v(X,Y,alpha,beta,gamma);
figure(1); s = surf(X,Y,Vactual); set(s,'LineStyle','none');
Vapprox = solvePoissonRect(F,a,b,c,d);
figure(2); s = surf(X,Y,Vapprox); set(s,'LineStyle','none');
max(abs(Vapprox(:)-Vactual(:)))

function true = v(x,y,a,b,g)
true = y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end

function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end