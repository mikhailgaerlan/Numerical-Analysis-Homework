function [vapp,X,Y,k] = alt_Schwarz(m,nmax,as,bs,tol,f,v)
%Multiplicative Schwarz Method
%   input  - m        grid length 1/2^m
%            nmax     max iterations
%            as       second rectangle x-axis starting point
%            bs       second rectangle y-axis starting point
%            tol      convergence tolerance
%            A        Laplacian
%            b        right-hand side vector
%            I1       matrix that returns region 1 elements
%            I2       matrix that returns region 2 elements
%   output - vapp     approximate solution
%            k        number of iterations reached
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

is_in_rectangle1 = rectangle1_x1 < X1 & X1 < rectangle1_x2 &...
    rectangle1_y1 < Y1 & Y1 < rectangle1_y2;
is_in_rectangle2 = rectangle2_x1 < X2 & X2 < rectangle2_x2 &...
    rectangle2_y1 < Y2 & Y2 < rectangle2_y2;
is_gamma1 = (X1 == rectangle1_x2 & Y1 > rectangle2_y1 & Y1 < rectangle2_y2);
is_gamma2 = (X2 == rectangle2_x1 & Y2 < rectangle2_y2 & Y2 > rectangle2_y1);
is_gamma3 = (Y2 == rectangle2_y1 & X2 < rectangle1_x2 & X2 > rectangle2_x1);
is_gamma4 = (Y2 == rectangle2_y2 & X2 < rectangle1_x2 & X2 > rectangle2_x1);
is_next_to_gamma1 = (X1 == rectangle1_x2-h & Y1 > rectangle2_y1 & Y1 < rectangle2_y2);
is_next_to_gamma2 = (X2 == rectangle2_x1+h & Y2 < rectangle2_y2 & Y2 > rectangle2_y1);
is_next_to_gamma3 = (Y2 == rectangle2_y1+h & X2 < rectangle1_x2 & X2 > rectangle2_x1);
is_next_to_gamma4 = (Y2 == rectangle2_y2-h & X2 < rectangle1_x2 & X2 > rectangle2_x1);

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
to_ignore2 = (X2 < rectangle2_x1+1/5);
to_accept1 = is_in_rectangle1;
to_accept2 = is_in_rectangle2 & ~to_ignore2;

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
ind_to_accept1 = get_indices(to_accept1);
ind_to_accept2 = get_indices(to_accept2);

v1 = reshape(v(X1,Y1),[A1_m,1]);
v1(ind_in_rectangle1) = zeros(length(ind_in_rectangle1),1);
v1(ind_gamma1) = zeros(length(ind_gamma1),1);
v2 = reshape(v(X2,Y2),[A2_m,1]);
v2(ind_in_rectangle2) = zeros(length(ind_in_rectangle2),1);
v2(ind_gamma2) = zeros(length(ind_gamma2),1);
v2(ind_gamma3) = zeros(length(ind_gamma3),1);
v2(ind_gamma4) = zeros(length(ind_gamma4),1);

b1 = reshape(h^2*f(X1,Y1),[A1_m,1]);
b1(ind_next_to_rect1_x1) = b1(ind_next_to_rect1_x1) + v1(ind_rect1_x1);
b1(ind_next_to_rect1_x2) = b1(ind_next_to_rect1_x2) + v1(ind_rect1_x2);
b1(ind_next_to_rect1_y1) = b1(ind_next_to_rect1_y1) + v1(ind_rect1_y1);
b1(ind_next_to_rect1_y2) = b1(ind_next_to_rect1_y2) + v1(ind_rect1_y2);
b2 = reshape(h^2*f(X2,Y2),[A2_m,1]);
b2(ind_next_to_rect2_x2) = b2(ind_next_to_rect2_x2) + v2(ind_rect2_x2);
b2(ind_next_to_rect2_y1) = b2(ind_next_to_rect2_y1) + v2(ind_rect2_y1);
b2(ind_next_to_rect2_y2) = b2(ind_next_to_rect2_y2) + v2(ind_rect2_y2);

for n = 1:nmax
    v1(ind_in_rectangle1) = A1inv(b1 + bgamma1Igamma1region2(v2));
    v2(ind_in_rectangle2) = A2inv(b2 + bgamma2Igamma2region1(v1));
    k = n;
end
vapp = [v1(ind_to_accept1);v2(ind_to_accept2)];
X = [X1(ind_to_accept1);X2(ind_to_accept2)];
Y = [Y1(ind_to_accept1);Y2(ind_to_accept2)];
    function b1 = bgamma1Igamma1region2(v2)
        b1 = zeros(A1_m,1); b1p = zeros(A1_m,1);
        for i = ind_gamma1'
            [x,y,z] = get_four_points(v2,i,X1_reshape,Y1_reshape,X2,Y2,h);
            b1p(i) = interpolate(X1_reshape(i),Y1_reshape(i),x,y,z);
        end; b1(ind_next_to_gamma1) = b1p(ind_gamma1);
    end
    function b2 = bgamma2Igamma2region1(v1)
        b2 = zeros(A2_m,1); b2p = zeros(A2_m,1);
        for i = ind_gamma2'
            [x,y,z] = get_four_points(v1,i,X2_reshape,Y2_reshape,X1,Y1,h);
            b2p(i) = interpolate(X2_reshape(i),Y2_reshape(i),x,y,z);
        end; b2(ind_next_to_gamma2) = b2p(ind_gamma2);
        for i = ind_gamma3'
            [x,y,z] = get_four_points(v1,i,X2_reshape,Y2_reshape,X1,Y1,h);
            b2p(i) = interpolate(X2_reshape(i),Y2_reshape(i),x,y,z);
        end; b2(ind_next_to_gamma3) = b2p(ind_gamma3);
        for i = ind_gamma4'
            [x,y,z] = get_four_points(v1,i,X2_reshape,Y2_reshape,X1,Y1,h);
            b2p(i) = interpolate(X2_reshape(i),Y2_reshape(i),x,y,z);
        end; b2(ind_next_to_gamma4) = b2p(ind_gamma4);
    end
    function b = A1inv(v)
        m1 = rectangle1_m1-2; m2 = rectangle1_m2-2;
        z = reshape(v(ind_in_rectangle1),[m2,m1]);
        g = solvePoisson(z,rectangle1_x1,rectangle1_x2,...
            rectangle1_y1,rectangle1_y2);
        b = reshape(g,[m1*m2,1]);
    end
    function b = A2inv(v)
        m1 = rectangle2_m1-2; m2 = rectangle2_m2-2;
        z = reshape(v(ind_in_rectangle2),[m2,m1]);
        g = solvePoisson(z,rectangle2_x1,rectangle2_x2,...
            rectangle2_y1,rectangle2_y2);
        b = reshape(g,[m1*m2,1]);
    end
    function [x,y,z] = get_four_points(v2,ind,X1,Y1,X2,Y2,h)
        is_near_point = X2 > X1(ind)-h & X2 < X1(ind)+h &...
            Y2 > Y1(ind)-h & Y2 < Y1(ind)+h;
        ind_near_point = get_indices(is_near_point);
        x = X2(ind_near_point); y = Y2(ind_near_point); z = v2(ind_near_point);
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
        nn = size(I,2);
        pre_ind = I'*I*(1:nn)';
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
end