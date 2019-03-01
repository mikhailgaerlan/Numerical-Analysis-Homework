clearvars; clc; m = 4; nmax = 5;
as = 2/5; bs = 9/5; alpha = 2; beta = 7/2; gamma = 2;

rectangle1_x1=0;rectangle1_x2=1; rectangle1_y1=0; rectangle1_y2=3;
rectangle2_x1=as;rectangle2_x2=as+4;rectangle2_y1=bs; rectangle2_y2=bs+1;
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

is_rect1_x1 = (X1 == rectangle1_x1 & Y1 > rectangle1_y1 & Y1 < rectangle1_y2);
is_rect1_x2 = (X1 == rectangle1_x2 & Y1 > rectangle1_y1 & Y1 < rectangle1_y2 &...
    (Y1 > rectangle2_y2 | Y1 < rectangle2_y1));
is_rect1_y1 = (Y1 == rectangle1_y1 & X1 > rectangle1_x1 & X1 < rectangle1_x2);
is_rect1_y2 = (Y1 == rectangle1_y2 & X1 > rectangle1_x1 & X1 < rectangle1_x2);
is_rect2_x2 = (X2 == rectangle2_x2 & Y2 > rectangle2_y1 & Y2 < rectangle2_y2);
is_rect2_y1 = (Y2 == rectangle2_y1 & X2 > rectangle2_x1 &...
    X2 > rectangle1_x2 & Y2 < rectangle2_x2);
is_rect2_y2 = (Y2 == rectangle2_y2 & X2 > rectangle2_x1 &...
    X2 > rectangle1_x2 & Y2 < rectangle2_x2);

is_next_to_rect1_x1 = (X1 == rectangle1_x1+h & Y1 > rectangle1_y1 & Y1 < rectangle1_y2);
is_next_to_rect1_x2 = (X1 == rectangle1_x2-h & Y1 > rectangle1_y1 & Y1 < rectangle1_y2 &...
    (Y1 > rectangle2_y2 | Y1 < rectangle2_y1));
is_next_to_rect1_y1 = (Y1 == rectangle1_y1+h & X1 > rectangle1_x1 & X1 < rectangle1_x2);
is_next_to_rect1_y2 = (Y1 == rectangle1_y2-h & X1 > rectangle1_x1 & X1 < rectangle1_x2);
is_next_to_rect2_x2 = (X2 == rectangle2_x2-h & Y2 > rectangle2_y1 & Y2 < rectangle2_y2);
is_next_to_rect2_y1 = (Y2 == rectangle2_y1+h & X2 > rectangle2_x1 &...
    X2 > rectangle1_x2 & Y2 < rectangle2_x2);
is_next_to_rect2_y2 = (Y2 == rectangle2_y2-h & X2 > rectangle2_x1 &...
    X2 > rectangle1_x2 & Y2 < rectangle2_x2);

ind_rect1_x1 = get_indices(is_rect1_x1);
ind_rect1_x2 = get_indices(is_rect1_x2);
ind_rect1_y1 = get_indices(is_rect1_y1);
ind_rect1_y2 = get_indices(is_rect1_y2);
ind_rect2_x2 = get_indices(is_rect2_x2);
ind_rect2_y1 = get_indices(is_rect2_y1);
ind_rect2_y2 = get_indices(is_rect2_y2);

ind_next_to_rect1_x1 = get_indices(is_next_to_rect1_x1);
ind_next_to_rect1_x2 = get_indices(is_next_to_rect1_x2);
ind_next_to_rect1_y1 = get_indices(is_next_to_rect1_y1);
ind_next_to_rect1_y2 = get_indices(is_next_to_rect1_y2);
ind_next_to_rect2_x2 = get_indices(is_next_to_rect2_x2);
ind_next_to_rect2_y1 = get_indices(is_next_to_rect2_y1);
ind_next_to_rect2_y2 = get_indices(is_next_to_rect2_y2);

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

plot(X1(is_in_rectangle1),Y1(is_in_rectangle1),'o'); hold on
plot(X2(is_in_rectangle2),Y2(is_in_rectangle2),'o');
plot(X2(is_rect2_x2),Y2(is_rect2_x2),'o');
plot(X2(is_next_to_rect2_x2),Y2(is_next_to_rect2_x2),'o'); hold off

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