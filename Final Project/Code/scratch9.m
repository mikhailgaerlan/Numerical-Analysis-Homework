clearvars; clc; m = 3; nmax = 5;
a = 2/5; b = 9/5; alpha = 0; beta = 1; gamma = 1/2;

rectangle1_x1 = 0; rectangle1_x2 = 1;
rectangle1_y1 = 0; rectangle1_y2 = 3;

rectangle2_x1 = a; rectangle2_x2 = a+4;
rectangle2_y1 = b; rectangle2_y2 = b+1;

rectangle1_m1 = 1*2^m+1; rectangle1_m2 = 3*2^m+1;
rectangle2_m1 = 4*2^m+1; rectangle2_m2 = 1*2^m+1;

h = 1/2^m;

[X1,Y1]=meshgrid(...
    rectangle1_x1:h:rectangle1_x2,...
    rectangle1_y1:h:rectangle1_y2);

[X2,Y2]=meshgrid(...
    rectangle2_x1:h:rectangle2_x2,...
    rectangle2_y1:h:rectangle2_y2);

is_rectangle1 = true(rectangle1_m2,rectangle1_m1);
is_rectangle2 = true(rectangle2_m2,rectangle2_m1);

r1_is_in_rectangle1 = rectangle1_x1 < X1 & X1 < rectangle1_x2 &...
    rectangle1_y1 < Y1 & Y1 < rectangle1_y2;
r2_is_in_rectangle2 = rectangle2_x1 < X2 & X2 < rectangle2_x2 &...
    rectangle2_y1 < Y2 & Y2 < rectangle2_y2;
r1_is_in_border1 = (...
    (rectangle1_x1 == X1 & rectangle1_y1 < Y1 & Y1 < rectangle1_y2) | ...
    (rectangle1_x2 == X1 & rectangle1_y1 < Y1 & Y1 < rectangle1_y2 &...
    (Y1 > rectangle2_y2 | Y1 < rectangle2_y1)) | ...
    (rectangle1_y1 == Y1 & rectangle1_x1 < X1 & X1 < rectangle1_x2) | ...
    (rectangle1_y2 == Y1 & rectangle1_x1 < X1 & X1 < rectangle1_x2));
r2_is_in_border2 = (...
    (rectangle2_x2 == X2 & rectangle2_y1 < Y2 & Y2 < rectangle2_y2 &...
    X2 > rectangle1_x2) | (rectangle2_y1 == Y2 & rectangle2_x1 < X2 &...
    X2 < rectangle2_x2 & X2 > rectangle1_x2) | (rectangle2_y2 == Y2 &...
    rectangle2_x1 < X2 & X2 < rectangle2_x2 & X2 > rectangle1_x2));
r1_is_in_border2 = (X1 == rectangle1_x2 & Y1 > rectangle2_y1 & Y1 < rectangle2_y2);
r2_is_in_border1 = (...
    (X2 == rectangle2_x1 & Y2 < rectangle2_y2 & Y2 > rectangle2_y1) |...
    (Y2 == rectangle2_y1 & X2 < rectangle1_x2 & X2 > rectangle2_x1) |...
    (Y2 == rectangle2_y2 & X2 < rectangle1_x2 & X2 > rectangle2_x1));

figure(1); plot(X1(r1_is_in_rectangle1),Y1(r1_is_in_rectangle1),'o'); hold on
plot(X2(r2_is_in_rectangle2),Y2(r2_is_in_rectangle2),'o');
plot(X1(r1_is_in_border1),Y1(r1_is_in_border1),'o');
plot(X2(r2_is_in_border2),Y2(r2_is_in_border2),'o');
plot(X1(r1_is_in_border2),Y1(r1_is_in_border2),'o');
plot(X2(r2_is_in_border1),Y2(r2_is_in_border1),'o'); hold off

figure(2);
plot(X1(r1_is_in_border1),Y1(r1_is_in_border1),'o'); hold on
plot(X2(r2_is_in_border2),Y2(r2_is_in_border2),'o');
plot(X1(r1_is_in_border2),Y1(r1_is_in_border2),'o');
plot(X2(r2_is_in_border1),Y2(r2_is_in_border1),'o'); hold off

r1_is_next_to_1_x2 = (X1 == rectangle1_x2-h &...
    Y1 < rectangle2_y2 & Y1 > rectangle2_y1);

figure(3);
plot(X1(r1_is_in_border1),Y1(r1_is_in_border1),'o'); hold on
plot(X2(r2_is_in_border2),Y2(r2_is_in_border2),'o');
plot(X1(r1_is_next_to_1_x2),Y1(r1_is_next_to_1_x2),'o');
plot(X2(r2_is_in_border2),Y2(r2_is_in_border2),'o');
plot(X1(r1_is_in_border2),Y1(r1_is_in_border2),'o'); hold off

r2_is_next_to_2_x1 = (...
    X2 == rectangle2_x1+h & Y2 > rectangle2_y1 & Y2 < rectangle2_y2);
r2_is_next_to_2_y1 = (...
    Y2 == rectangle2_y1+h & X2 > rectangle2_x1 & X2 < rectangle1_x2);
r2_is_next_to_2_y2 = (...
    Y2 == rectangle2_y2-h & X2 > rectangle2_x1 & X2 < rectangle1_x2);

figure(4);
plot(X1(r1_is_in_border1),Y1(r1_is_in_border1),'o'); hold on
plot(X2(r2_is_in_border2),Y2(r2_is_in_border2),'o');
plot(X2(r2_is_next_to_2_x1),Y2(r2_is_next_to_2_x1),'o');
plot(X2(r2_is_next_to_2_y1),Y2(r2_is_next_to_2_y1),'o');
plot(X2(r2_is_next_to_2_y2),Y2(r2_is_next_to_2_y2),'o');
plot(X2(r2_is_in_border2),Y2(r2_is_in_border2),'o');
plot(X1(r1_is_in_border2),Y1(r1_is_in_border2),'o'); hold off

%r2_is_interp_x2 = (...
%    Y2 == rectangle2_y2-h & X2 > rectangle2_x1 & X2 < rectangle1_x2);

%is_next_to_x2 = X==rectangle1_x2-h & rectangle1_y1<Y &...
%    Y<rectangle1_y2 & ~is_in_rectangle1_and_2;
%is_next_to_y1 = Y==rectangle1_y1+h & rectangle1_x1<X & X<rectangle1_x2;
%is_next_to_y2 = Y==rectangle1_y2-h & rectangle1_x1<X & X<rectangle1_x2;
%is_next_to_b = Y==rectangle2_y1+h & rectangle2_x1<X &...
%    X<rectangle2_x2 & ~is_in_rectangle1_and_2 & is_in_rectangle2;
%is_next_to_b_1 = Y==rectangle2_y2-h & rectangle2_x1<X &...
%    X<rectangle2_x2 & ~is_in_rectangle1_and_2 & is_in_rectangle2;
%is_next_to_a_4 = X==rectangle2_x2-h & rectangle2_y1<Y & Y<rectangle2_y2;

%I1_get_border_1_pre = is_in_rectangle1_reshape(is_in_rectangle12_reshape);
%I1_elems = sum(I1_pre); otherI1 = I1_pre.*(1:total_elems)';
%ind1 = otherI1(otherI1 ~= 0);
%I1 = sparse(1:I1_elems,ind1,ones(I1_elems,1),I1_elems,total_elems);

I1_rectangle1 = I_region(r1_is_in_rectangle1);
figure(5);spy(I1_rectangle1'*I1_rectangle1)

I1_border1 = I_region(r1_is_in_border1);
figure(6);spy(I1_border1'*I1_border1)

F1 = f(X1,Y1,alpha,beta,gamma); G1 = v(X1,Y1,alpha,beta,gamma);
H1 = zeros(rectangle1_m2,rectangle1_m1);
reshape(G1,[rectangle1_m1*rectangle1_m2,1]);

I2_border2 = I_region(r2_is_in_border2);
figure(6);spy(I2_border2'*I2_border2);

F2 = f(X2,Y2,alpha,beta,gamma); G2 = v(X2,Y2,alpha,beta,gamma);
H2 = zeros(rectangle2_m2,rectangle2_m1);
indices = get_indices(r1_is_in_border1);
G2_reshape = reshape(G2,[rectangle2_m1*rectangle2_m2,1]);
G2_reshape(indices)

%A1_m = rectangle1_m1*rectangle1_m2;
%r1_is_in_rectangle1_reshape = reshape(r1_is_in_rectangle1,[rectangle1_m1*rectangle1_m2,1]);
%is_in_rectangle1_reshape = reshape(r1_is_in_rectangle1,[rectangle1_m1*rectangle1_m2,1]);
%total_elems1 = sum(is_in_rectangle1_reshape);
%I1_in_rectangle1_pre = r1_is_in_rectangle1_reshape(is_rectangle1);
%I1_elems = sum(I1_in_rectangle1_pre);
%otherI1 = I1_in_rectangle1_pre.*(1:total_elems1)';
%ind1 = otherI1(otherI1 ~= 0);
%I1_in_rectangle1 = sparse(1:I1_elems,ind1,ones(I1_elems,1),I1_elems,total_elems1);
%I1_in_rectangle1*(1:A1_m')
%figure(5);spy(I1_in_rectangle1)

%Create the matrices that get the elements of regions
%total_elems = sum(is_in_rectangle12_reshape);
%I1_pre = is_in_rectangle1_reshape(is_in_rectangle12_reshape);
%I1_elems = sum(I1_pre); otherI1 = I1_pre.*(1:total_elems)';
%ind1 = otherI1(otherI1 ~= 0);
%I1 = sparse(1:I1_elems,ind1,ones(I1_elems,1),I1_elems,total_elems);
%I2_pre = is_in_rectangle2_reshape(is_in_rectangle12_reshape);
%I2_elems = sum(I2_pre); otherI2 = I2_pre.*(1:total_elems)';
%ind2 = otherI2(otherI2 ~= 0);
%I2 = sparse(1:I2_elems,ind2,ones(I2_elems,1),I2_elems,total_elems);

%Create the right hand-side b1
%F1 = f(X1,Y1,alpha,beta,gamma); G1 = v(X1,Y1,alpha,beta,gamma);
%H = zeros(rectangle1_m2,rectangle1_m1);
%for i = 1:m1
%    for j = 1:m2
%        if is_in_rectangle1_or_2(j,i); H(j,i) = h^2*F(j,i); end
%        if is_next_to_x1(j,i); H(j,i) = H(j,i)  + G(j,i-1); end
%        if is_next_to_x2(j,i); H(j,i) = H(j,i)  + G(j,i+1); end
%        if is_next_to_y1(j,i); H(j,i) = H(j,i)  + G(j-1,i); end
%        if is_next_to_y2(j,i); H(j,i) = H(j,i)  + G(j+1,i); end
%        if is_next_to_b(j,i); H(j,i) = H(j,i)   + G(j-1,i); end
%        if is_next_to_b_1(j,i); H(j,i) = H(j,i) + G(j+1,i); end
%        if is_next_to_a_4(j,i); H(j,i) = H(j,i) + G(j,i+1); end
%    end
%end
%b = H(is_in_rectangle1_or_2);
%X1s = X(is_in_rectangle1_or_2);Y1s = Y(is_in_rectangle1_or_2);

function indices = get_indices(logical_region)
I = I_region(logical_region);
n = size(I,2);
pre_ind = I'*I*(1:n)';
indices = pre_ind(pre_ind~=0);
end

function I = I_region(logical_region)
[m1,m2] = size(logical_region);
A_m = m1*m2;
logical_reshape = reshape(logical_region,[A_m,1]);
num_elems = sum(logical_reshape);
ind1 = logical_reshape.*(1:A_m)';
ind = ind1(ind1~=0);
I = sparse(1:num_elems,ind,ones(num_elems,1),num_elems,A_m,num_elems);
end
function true = v(x,y,a,b,g); true=y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y); end
function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end