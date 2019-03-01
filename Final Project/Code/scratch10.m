clearvars; clc; m = 4; nmax = 5;
a = 2/5; b = 9/5; alpha = 2; beta = 7/2; gamma = 2;

rectangle1_x1 = 0; rectangle1_x2 = 1;
rectangle1_y1 = 0; rectangle1_y2 = 3;
rectangle2_x1 = a; rectangle2_x2 = a+4;
rectangle2_y1 = b; rectangle2_y2 = b+1;
rectangle1_m1 = 1*2^m+1; rectangle1_m2 = 3*2^m+1;
rectangle2_m1 = 4*2^m+1; rectangle2_m2 = 1*2^m+1;
A1_m = rectangle1_m1*rectangle1_m2; A2_m = rectangle2_m1*rectangle2_m2;

h = 1/2^m;

[X1,Y1]=meshgrid(...
    rectangle1_x1:h:rectangle1_x2,...
    rectangle1_y1:h:rectangle1_y2);
[X2,Y2]=meshgrid(...
    rectangle2_x1:h:rectangle2_x2,...
    rectangle2_y1:h:rectangle2_y2);
X1_reshape = reshape(X1,[A1_m,1]); Y1_reshape = reshape(Y1,[A1_m,1]);
X2_reshape = reshape(X2,[A2_m,1]); Y2_reshape = reshape(Y2,[A2_m,1]);

is_rectangle1 = true(rectangle1_m2,rectangle1_m1);
is_rectangle2 = true(rectangle2_m2,rectangle2_m1);
is_in_rectangle1 = rectangle1_x1 < X1 & X1 < rectangle1_x2 &...
    rectangle1_y1 < Y1 & Y1 < rectangle1_y2;
is_in_rectangle2 = rectangle2_x1 < X2 & X2 < rectangle2_x2 &...
    rectangle2_y1 < Y2 & Y2 < rectangle2_y2;
is_sigma1 = (...
    (rectangle1_x1 == X1 & rectangle1_y1 < Y1 & Y1 < rectangle1_y2) | ...
    (rectangle1_x2 == X1 & rectangle1_y1 < Y1 & Y1 < rectangle1_y2 &...
    (Y1 > rectangle2_y2 | Y1 < rectangle2_y1)) | ...
    (rectangle1_y1 == Y1 & rectangle1_x1 < X1 & X1 < rectangle1_x2) | ...
    (rectangle1_y2 == Y1 & rectangle1_x1 < X1 & X1 < rectangle1_x2));
is_sigma2 = (...
    (rectangle2_x2 == X2 & rectangle2_y1 < Y2 & Y2 < rectangle2_y2 &...
    X2 > rectangle1_x2) | (rectangle2_y1 == Y2 & rectangle2_x1 < X2 &...
    X2 < rectangle2_x2 & X2 > rectangle1_x2) | (rectangle2_y2 == Y2 &...
    rectangle2_x1 < X2 & X2 < rectangle2_x2 & X2 > rectangle1_x2));
is_gamma1 = (...
    X1 == rectangle1_x2 & Y1 > rectangle2_y1 & Y1 < rectangle2_y2);
is_gamma2 = (...
    X2 == rectangle2_x1 & Y2 < rectangle2_y2 & Y2 > rectangle2_y1);
is_gamma3 = (...
    Y2 == rectangle2_y1 & X2 < rectangle1_x2 & X2 > rectangle2_x1);
is_gamma4 = (...
    Y2 == rectangle2_y2 & X2 < rectangle1_x2 & X2 > rectangle2_x1);
is_next_to_gamma1 = (...
    X1 == rectangle1_x2-h & Y1 > rectangle2_y1 & Y1 < rectangle2_y2);
is_next_to_gamma2 = (...
    X2 == rectangle2_x1+h & Y2 < rectangle2_y2 & Y2 > rectangle2_y1);
is_next_to_gamma3 = (...
    Y2 == rectangle2_y1+h & X2 < rectangle1_x2 & X2 > rectangle2_x1);
is_next_to_gamma4 = (...
    Y2 == rectangle2_y2-h & X2 < rectangle1_x2 & X2 > rectangle2_x1);

ind_in_rectangle1 = get_indices(is_in_rectangle1);
ind_gamma1 = get_indices(is_gamma1);
ind_next_to_gamma1 = get_indices(is_next_to_gamma1);
ind_in_rectangle2 = get_indices(is_in_rectangle2);
ind_gamma2 = get_indices(is_gamma2);
ind_gamma3 = get_indices(is_gamma3);
ind_gamma4 = get_indices(is_gamma4);
ind_next_to_gamma2 = get_indices(is_next_to_gamma2);
ind_next_to_gamma3 = get_indices(is_next_to_gamma3);
ind_next_to_gamma4 = get_indices(is_next_to_gamma4);

v1 = reshape(v(X1,Y1,alpha,beta,gamma),[A1_m,1]);
v1(ind_in_rectangle1) = zeros(length(ind_in_rectangle1),1);
v1(ind_gamma1) = zeros(length(ind_gamma1),1);

v2 = reshape(v(X2,Y2,alpha,beta,gamma),[A2_m,1]);
%v2(ind_in_rectangle2) = zeros(length(ind_in_rectangle2),1);
v2(ind_gamma2) = zeros(length(ind_gamma2),1);
v2(ind_gamma3) = zeros(length(ind_gamma3),1);
v2(ind_gamma4) = zeros(length(ind_gamma4),1);

figure(1); scatter3(X1_reshape,Y1_reshape,v1,'o','CData',v1); hold on
scatter3(X2_reshape,Y2_reshape,v2,'o','CData',v2); view(2); hold off

for i = ind_gamma1'
    [x,y,z] = get_four_points(v2,i,X1_reshape,Y1_reshape,X2,Y2,h);
    v1(i) = interpolate(X1_reshape(i),Y1_reshape(i),x,y,z);
end
v1(ind_next_to_gamma1) = v1(ind_next_to_gamma1) - v1(ind_gamma1);

figure(2); scatter3(X1_reshape,Y1_reshape,v1,'o','CData',v1); hold on
scatter3(X2_reshape,Y2_reshape,v2,'o','CData',v2); view(2);hold off

function [x,y,z] = get_four_points(v2,ind,X1,Y1,X2,Y2,h)
is_near_point = X2 > X1(ind)-h & X2 < X1(ind)+h &...
    Y2 > Y1(ind)-h & Y2 < Y1(ind)+h;
ind_near_point = get_indices(is_near_point);
x = X2(ind_near_point); y = X2(ind_near_point); z = v2(ind_near_point);
end

function u = interpolate(x0,y0,x,y,z)
interpmatrix = [1,x(1),y(1),x(1)*y(1);...
    1,x(2),y(2),x(2)*y(2);...
    1,x(3),y(3),x(3)*y(3);...
    1,x(4),y(4),x(4)*y(4)];
c = interpmatrix\z;
u = c(1)+c(2)*x0+c(3)*y0+c(4)*x0*y0;
end

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