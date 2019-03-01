clearvars; clc; m = 3; nmax = 5;
a = 1/2; b = 2; alpha = 2; beta = 7/2; gamma = 2;

rectangle1_x1 = 0; rectangle1_x2 = 1;
rectangle1_y1 = 0; rectangle1_y2 = 3;

rectangle2_x1 = a; rectangle2_x2 = a+4;
rectangle2_y1 = b; rectangle2_y2 = b+1;

h = 1/2^m;

%Beginning of function
rectangle_x1 = min([rectangle1_x1,rectangle2_x1]);
rectangle_x2 = max([rectangle1_x2,rectangle2_x2]);
rectangle_y1 = min([rectangle1_y1,rectangle2_y1]);
rectangle_y2 = max([rectangle1_y2,rectangle2_y2]);

[X,Y]=meshgrid(...
    rectangle_x1:h:rectangle_x2,...
    rectangle_y1:h:rectangle_y2);

m1 = (rectangle_x2-rectangle_x1)/h+1;
m2 = (rectangle_y2-rectangle_y1)/h+1;

is_in_rectangle1 = rectangle1_x1 < X & X < rectangle1_x2 &...
    rectangle1_y1 < Y & Y < rectangle1_y2;
is_in_rectangle2 = rectangle2_x1 < X & X < rectangle2_x2 &...
    rectangle2_y1 < Y & Y < rectangle2_y2;
is_in_rectangle1_or_2 = is_in_rectangle1 | is_in_rectangle2;
is_in_rectangle1_and_2 = is_in_rectangle1 & is_in_rectangle2;
%is_in_border1 = (...
%    (rectangle1_x1 == X & rectangle1_y1 < Y & Y < rectangle1_y2) | ...
%    (rectangle1_x2 == X & rectangle1_y1 < Y & Y < rectangle1_y2) | ...
%    (rectangle1_y1 == Y & rectangle1_x1 < X & X < rectangle1_x2) | ...
%    (rectangle1_y2 == Y & rectangle1_x1 < X & X < rectangle1_x2));
%is_in_border2 = (...
%    (rectangle2_x1 == X & rectangle2_y1 < Y & Y < rectangle2_y2) | ...
%    (rectangle2_x2 == X & rectangle2_y1 < Y & Y < rectangle2_y2) | ...
%    (rectangle2_y1 == Y & rectangle2_x1 < X & X < rectangle2_x2) | ...
%    (rectangle2_y2 == Y & rectangle2_x1 < X & X < rectangle2_x2));
%is_in_border = (is_in_border1 | is_in_border2) & (~is_in_rectangle1_or_2);

%figure(1)
%plot(X(is_in_rectangle1_or_2),Y(is_in_rectangle1_or_2),'o'); hold on
%plot(X(is_in_rectangle1_and_2),Y(is_in_rectangle1_and_2),'o'); hold on
%plot(X(is_in_border1),Y(is_in_border1),'o'); hold on
%plot(X(is_in_border2),Y(is_in_border2),'o'); hold off

%T_m1 = 2*speye(rectangle_m2)+...
%    sparse(2:rectangle_m2,1:rectangle_m2-1,-ones(rectangle_m2-1,1),rectangle_m2,rectangle_m2)+...
%    sparse(1:rectangle_m2-1,2:rectangle_m2,-ones(rectangle_m2-1,1),rectangle_m2,rectangle_m2);
%I_m1 = speye(rectangle_m2);
%A_m = rectangle_m1*rectangle_m2;
%A_rectangle = [T_m1+2*I_m1,-I_m1];
%for i=3:rectangle_m1; A_rectangle = [A_rectangle,z_m1];end
%rectangle_m1-2
%for i=2:(rectangle_m1-1)
%    tic
%    B = []; for j=1:(i-2); B = [B,z_m1]; end
%    B = [B,-I_m1,T_m1+2*I_m1,-I_m1];
%    for j=(i+2):rectangle_m1; B = [B,z_m1]; end
%    A_rectangle = [A_rectangle;B];
%    toc
%end
%B = []; for i=1:(rectangle_m1-2); B = [B,z_m1]; end; B = [B,-I_m1,T_m1+2*I_m1];
%A_rectangle = [A_rectangle;B];

T_m1 = 2*speye(m2)+...
    sparse(2:m2,1:m2-1,-ones(m2-1,1),m2,m2)+...
    sparse(1:m2-1,2:m2,-ones(m2-1,1),m2,m2);
I_m1 = speye(m2);
A_m = m1*m2;
A_rectangle = sparse(A_m,A_m);
A_rectangle(1:m2,1:m2) = T_m1+2*I_m1;
A_rectangle(1:m2,m2+1:2*m2) = -I_m1;
for i=2:m1-1
    A_rectangle((i-1)*m2+1:i*m2,(i-1)*m2+1:i*m2) = T_m1+2*I_m1;
    A_rectangle((i-1)*m2+1:i*m2,i*m2+1:(i+1)*m2) = -I_m1;
    A_rectangle((i-1)*m2+1:i*m2,(i-2)*m2+1:(i-1)*m2) = -I_m1;
end
A_rectangle((m1-1)*m2+1:m1*m2,(m1-1)*m2+1:m1*m2) = T_m1+2*I_m1;
A_rectangle((m1-1)*m2+1:m1*m2,(m1-2)*m2+1:(m1-1)*m2) = -I_m1;

%figure(2)
%plot(X(is_in_rectangle1_or_2),Y(is_in_rectangle1_or_2),'o'); hold on
%plot(X(is_in_rectangle1),Y(is_in_rectangle1),'o'); hold on
%plot(X(is_in_border),Y(is_in_border),'o'); hold off

is_in_rectangle1_reshape = reshape(is_in_rectangle1,[A_m,1]);
is_in_rectangle2_reshape = reshape(is_in_rectangle2,[A_m,1]);
is_in_rectangle12_reshape = reshape(is_in_rectangle1_or_2,[A_m,1]);
A1 = A_rectangle(is_in_rectangle1_reshape,is_in_rectangle1_reshape);
A2 = A_rectangle(is_in_rectangle2_reshape,is_in_rectangle2_reshape);
A = A_rectangle(is_in_rectangle12_reshape,is_in_rectangle12_reshape);

total_elems = sum(is_in_rectangle12_reshape);
I1_pre = is_in_rectangle1_reshape(is_in_rectangle12_reshape);
I1_elems = sum(I1_pre);
otherI1 = I1_pre.*(1:total_elems)';
ind1 = otherI1(otherI1 ~= 0);
I1 = sparse(1:I1_elems,ind1,ones(I1_elems,1),I1_elems,total_elems);
I2_pre = is_in_rectangle2_reshape(is_in_rectangle12_reshape);
I2_elems = sum(I2_pre);
otherI2 = I2_pre.*(1:total_elems)';
ind2 = otherI2(otherI2 ~= 0);
I2 = sparse(1:I2_elems,ind2,ones(I2_elems,1),I2_elems,total_elems);

is_next_to_x1 = X==rectangle1_x1+h & rectangle1_y1<Y & Y<rectangle1_y2;
is_next_to_x2 = X==rectangle1_x2-h & rectangle1_y1<Y & Y<rectangle1_y2 & ~is_in_rectangle1_and_2;
is_next_to_y1 = Y==rectangle1_y1+h & rectangle1_x1<X & X<rectangle1_x2;
is_next_to_y2 = Y==rectangle1_y2-h & rectangle1_x1<X & X<rectangle1_x2;
is_next_to_b = Y==rectangle2_y1+h & rectangle2_x1<X & X<rectangle2_x2 & ~is_in_rectangle1_and_2 & is_in_rectangle2;
is_next_to_b_1 = Y==rectangle2_y2-h & rectangle2_x1<X & X<rectangle2_x2 & ~is_in_rectangle1_and_2 & is_in_rectangle2;
is_next_to_a_4 = X==rectangle2_x2-h & rectangle2_y1<Y & Y<rectangle2_y2;
%figure(3);
%plot(X(is_in_border),Y(is_in_border),'o'); hold on
%plot(X(is_next_to_x1),Y(is_next_to_x1),'o'); hold on
%plot(X(is_next_to_x2),Y(is_next_to_x2),'o'); hold on
%plot(X(is_next_to_y1),Y(is_next_to_y1),'o'); hold on
%plot(X(is_next_to_y2),Y(is_next_to_y2),'o'); hold on
%plot(X(is_next_to_b),Y(is_next_to_b),'o'); hold on
%plot(X(is_next_to_b_1),Y(is_next_to_b_1),'o'); hold on
%plot(X(is_next_to_a_4),Y(is_next_to_a_4),'o'); hold off


F = f(X,Y,alpha,beta,gamma);
G = v(X,Y,alpha,beta,gamma);
H = zeros(m2,m1);
for i = 1:m1
    for j = 1:m2
        if is_in_rectangle1_or_2(j,i); H(j,i) = h^2*F(j,i); end
        if is_next_to_x1(j,i); H(j,i) = H(j,i)  + G(j,i-1); end
        if is_next_to_x2(j,i); H(j,i) = H(j,i)  + G(j,i+1); end
        if is_next_to_y1(j,i); H(j,i) = H(j,i)  + G(j-1,i); end
        if is_next_to_y2(j,i); H(j,i) = H(j,i)  + G(j+1,i); end
        if is_next_to_b(j,i); H(j,i) = H(j,i)   + G(j-1,i); end
        if is_next_to_b_1(j,i); H(j,i) = H(j,i) + G(j+1,i); end
        if is_next_to_a_4(j,i); H(j,i) = H(j,i) + G(j,i+1); end
    end
end
%figure(5);spy(H);
%figure(6);scatter3(X(is_in_rectangle1_or_2),Y(is_in_rectangle1_or_2),H(is_in_rectangle1_or_2));
b = H(is_in_rectangle1_or_2);

vapprox = zeros(length(b),1);
vactual = v(X,Y,alpha,beta,gamma);
vactual = vactual(is_in_rectangle1_or_2);
%A1i = inv(A1);
%A2i = inv(A2);
for n = 1:nmax
    tic
    
    m1 = rectangle1_m1-2; m2 = rectangle1_m2-2;
    size(I1)
    size(I2)
    [m1,m2]
    z = reshape(I1*(b-A*vapprox),[m2,m1]);
    p = solvePoisson(z,rectangle1_x1,rectangle1_x2,rectangle1_y1,rectangle1_y2);
    vapprox = vapprox + I1'*reshape(p,[m1*m2,1]);
    
    m1 = rectangle2_m1-2; m2 = rectangle2_m2-2;
    z = reshape(I2*(b-A*vapprox),[m2,m1]);
    p = solvePoisson(z,rectangle2_x1,rectangle2_x2,rectangle2_y1,rectangle2_y2);
    vapprox = vapprox + I2'*reshape(p,[m1*m2,1]);
    %vapprox = vapprox + I1'*A1i*I1*(b-A*vapprox);
    %vapprox = vapprox + I2'*A2i*I2*(b-A*vapprox);
    toc
end
disp(max(abs(vapprox-vactual))/max(abs(vactual)))
figure(1);scatter3(X(is_in_rectangle1_or_2),Y(is_in_rectangle1_or_2),vactual,'square','CData',vactual);view(2);
figure(2);scatter3(X(is_in_rectangle1_or_2),Y(is_in_rectangle1_or_2),vapprox,'square','CData',vapprox);view(2);

function true = v(x,y,a,b,g)
true = y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end

function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end