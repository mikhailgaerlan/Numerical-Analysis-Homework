clearvars; clc; m = 7; nmax = 100;
a = 1/2; b = 2; alpha = 2; beta = 7/2; gamma = 2;

rectangle1_x1 = 0; rectangle1_x2 = 1;
rectangle1_y1 = 0; rectangle1_y2 = 3;

rectangle2_x1 = a; rectangle2_x2 = a+4;
rectangle2_y1 = b; rectangle2_y2 = b+1;

rectangle1_m1 = 1*2^m+1; rectangle1_m2 = 3*2^m+1;
rectangle2_m1 = 4*2^m+1; rectangle2_m2 = 1*2^m+1;

tic
[A,A1,A2,I1,I2,X,Y,b] = laplace_matching(m,a,b,@(x,y) f(x,y,alpha,beta,gamma),@(x,y) v(x,y,alpha,beta,gamma));
toc

vapprox = zeros(length(b),1);
vactual = v(X,Y,alpha,beta,gamma);
error = b - A*vapprox;
tic
for n = 1:nmax
    p = multB(b-A*vapprox,I1,rectangle1_m1,rectangle1_m2,...
        rectangle1_x1,rectangle1_x2,rectangle1_y1,rectangle1_y2);
    vapprox = vapprox + p;
    
    p = multB(b-A*vapprox,I2,rectangle2_m1,rectangle2_m2,...
        rectangle2_x1,rectangle2_x2,rectangle2_y1,rectangle2_y2);
    vapprox = vapprox + p;
    
    error = error-Minv(A*error,A,I1,rectangle1_m1,rectangle1_m2,...
        rectangle1_x1,rectangle1_x2,rectangle1_y1,rectangle1_y2,...
        I2,rectangle2_m1,rectangle2_m2,...
        rectangle2_x1,rectangle2_x2,rectangle2_y1,rectangle2_y2);
    if norm(error) < 1e-15; disp(n); break; end
end
toc
disp(max(abs(vapprox-vactual))/max(abs(vactual)))
figure(1);scatter3(X,Y,vactual,'square','CData',vactual);view(2);
figure(2);scatter3(X,Y,vapprox,'square','CData',vapprox);view(2);

function b = multB(v,I1,m1,m2,x1,x2,y1,y2)
m1p = m1-2; m2p = m2-2;
z = reshape(I1*v,[m2p,m1p]);
p = solvePoisson(z,x1,x2,y1,y2);
b = I1'*reshape(p,[m1p*m2p,1]);
end

function b=Minv(v,A,I1,m11,m12,x11,x12,y11,y12,I2,m21,m22,x21,x22,y21,y22)
b = multB(v,I1,m11,m12,x11,x12,y11,y12)+...
    multB(v,I2,m21,m22,x21,x22,y21,y22)-...
    multB(A*multB(v,I1,m11,m12,x11,x12,y11,y12),...
    I2,m21,m22,x21,x22,y21,y22);
end

function true = v(x,y,a,b,g)
true = y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end

function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end