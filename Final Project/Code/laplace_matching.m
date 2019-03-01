function [A,A1,A2,I1,I2,Xs,Ys,b] = laplace_matching(m,as,bs,f,v)
%Create Laplacian for two rectangles with matching grids Av=b
%   input  - as       second rectangle x-axis starting point
%            bs       second rectangle y-axis starting point
%            m        grid length 1/2^m
%   output - A        laplacian
%            A1       region 1 laplacian
%            A2       region 2 laplacian
%            I1       matrix that returns region 1 elements
%            I2       matrix that returns region 2 elements
%            Xs       x-values of regions 1 and 2
%            Ys       y-values of regions 1 and 2
%            b        b such that Av = b

%Parameter initialization
h = 1/2^m;
rectangle1_x1=0;rectangle1_x2=1; rectangle1_y1=0; rectangle1_y2=3;
rectangle2_x1=as;rectangle2_x2=as+4;rectangle2_y1=bs; rectangle2_y2=bs+1;
rectangle_x1 = min([rectangle1_x1,rectangle2_x1]);
rectangle_x2 = max([rectangle1_x2,rectangle2_x2]);
rectangle_y1 = min([rectangle1_y1,rectangle2_y1]);
rectangle_y2 = max([rectangle1_y2,rectangle2_y2]);
m1 = (rectangle_x2-rectangle_x1)/h+1;
m2 = (rectangle_y2-rectangle_y1)/h+1;

%Create a large rectangular grid that includes both regions
[X,Y]=meshgrid(...
    rectangle_x1:h:rectangle_x2,...
    rectangle_y1:h:rectangle_y2);

%Create the Laplacian matrix for the large rectangular area
T_m1 = 2*speye(m2)+...
    sparse(2:m2,1:m2-1,-ones(m2-1,1),m2,m2)+...
    sparse(1:m2-1,2:m2,-ones(m2-1,1),m2,m2); A_m = m1*m2;
I_m1 = speye(m2); A_rectangle=spalloc(A_m,A_m,5*A_m-2*m1-2*m2);
A_rectangle(1:m2,1:m2) = T_m1+2*I_m1; A_rectangle(1:m2,m2+1:2*m2) = -I_m1;
for i=2:m1-1
    A_rectangle((i-1)*m2+1:i*m2,(i-1)*m2+1:i*m2) = T_m1+2*I_m1;
    A_rectangle((i-1)*m2+1:i*m2,i*m2+1:(i+1)*m2) = -I_m1;
    A_rectangle((i-1)*m2+1:i*m2,(i-2)*m2+1:(i-1)*m2) = -I_m1;
end
A_rectangle((m1-1)*m2+1:m1*m2,(m1-1)*m2+1:m1*m2) = T_m1+2*I_m1;
A_rectangle((m1-1)*m2+1:m1*m2,(m1-2)*m2+1:(m1-1)*m2) = -I_m1;

%Logicals that specify regions of the rectangular area
is_in_rectangle1 = rectangle1_x1 < X & X < rectangle1_x2 &...
    rectangle1_y1 < Y & Y < rectangle1_y2;
is_in_rectangle2 = rectangle2_x1 < X & X < rectangle2_x2 &...
    rectangle2_y1 < Y & Y < rectangle2_y2;
is_in_rectangle1_or_2 = is_in_rectangle1 | is_in_rectangle2;
is_in_rectangle1_and_2 = is_in_rectangle1 & is_in_rectangle2;
is_in_rectangle1_reshape = reshape(is_in_rectangle1,[A_m,1]);
is_in_rectangle2_reshape = reshape(is_in_rectangle2,[A_m,1]);
is_in_rectangle12_reshape = reshape(is_in_rectangle1_or_2,[A_m,1]);

%Filter out the regions of interest
A1 = A_rectangle(is_in_rectangle1_reshape,is_in_rectangle1_reshape);
A2 = A_rectangle(is_in_rectangle2_reshape,is_in_rectangle2_reshape);
A = A_rectangle(is_in_rectangle12_reshape,is_in_rectangle12_reshape);

%Create the matrices that get the elements of regions
total_elems = sum(is_in_rectangle12_reshape);
I1_pre = is_in_rectangle1_reshape(is_in_rectangle12_reshape);
I1_elems = sum(I1_pre); otherI1 = I1_pre.*(1:total_elems)';
ind1 = otherI1(otherI1 ~= 0);
I1 = sparse(1:I1_elems,ind1,ones(I1_elems,1),I1_elems,total_elems);
I2_pre = is_in_rectangle2_reshape(is_in_rectangle12_reshape);
I2_elems = sum(I2_pre); otherI2 = I2_pre.*(1:total_elems)';
ind2 = otherI2(otherI2 ~= 0);
I2 = sparse(1:I2_elems,ind2,ones(I2_elems,1),I2_elems,total_elems);

%Create the border logicals
is_next_to_x1 = X==rectangle1_x1+h & rectangle1_y1<Y & Y<rectangle1_y2;
is_next_to_x2 = X==rectangle1_x2-h & rectangle1_y1<Y &...
    Y<rectangle1_y2 & ~is_in_rectangle1_and_2;
is_next_to_y1 = Y==rectangle1_y1+h & rectangle1_x1<X & X<rectangle1_x2;
is_next_to_y2 = Y==rectangle1_y2-h & rectangle1_x1<X & X<rectangle1_x2;
is_next_to_b = Y==rectangle2_y1+h & rectangle2_x1<X &...
    X<rectangle2_x2 & ~is_in_rectangle1_and_2 & is_in_rectangle2;
is_next_to_b_1 = Y==rectangle2_y2-h & rectangle2_x1<X &...
    X<rectangle2_x2 & ~is_in_rectangle1_and_2 & is_in_rectangle2;
is_next_to_a_4 = X==rectangle2_x2-h & rectangle2_y1<Y & Y<rectangle2_y2;

%Create the right hand-side
F = f(X,Y); G = v(X,Y); H = zeros(m2,m1);
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
b = H(is_in_rectangle1_or_2);
Xs = X(is_in_rectangle1_or_2);Ys = Y(is_in_rectangle1_or_2);
end