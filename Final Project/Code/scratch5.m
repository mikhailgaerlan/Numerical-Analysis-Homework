clearvars; clc; m = 7; nmax = 5;
a = 1/4; b = 1; alpha = 2; beta = 7/2; gamma = 2;

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

rectangle1_m1 = (rectangle1_x2-rectangle1_x1)/h+1;
rectangle1_m2 = (rectangle1_y2-rectangle1_y1)/h+1;
rectangle2_m1 = (rectangle2_x2-rectangle2_x1)/h+1;
rectangle2_m2 = (rectangle2_y2-rectangle2_y1)/h+1;
m1 = (rectangle_x2-rectangle_x1)/h+1;
m2 = (rectangle_y2-rectangle_y1)/h+1;

is_in_rectangle1 = rectangle1_x1 < X & X < rectangle1_x2 &...
    rectangle1_y1 < Y & Y < rectangle1_y2;
is_in_rectangle2 = rectangle2_x1 < X & X < rectangle2_x2 &...
    rectangle2_y1 < Y & Y < rectangle2_y2;
is_in_rectangle1_or_2 = is_in_rectangle1 | is_in_rectangle2;
is_in_rectangle1_and_2 = is_in_rectangle1 & is_in_rectangle2;
is_in_border1 = (...
    (rectangle1_x1 == X & rectangle1_y1 < Y & Y < rectangle1_y2) | ...
    (rectangle1_x2 == X & rectangle1_y1 < Y & Y < rectangle1_y2) | ...
    (rectangle1_y1 == Y & rectangle1_x1 < X & X < rectangle1_x2) | ...
    (rectangle1_y2 == Y & rectangle1_x1 < X & X < rectangle1_x2));
is_in_border2 = (...
    (rectangle2_x1 == X & rectangle2_y1 < Y & Y < rectangle2_y2) | ...
    (rectangle2_x2 == X & rectangle2_y1 < Y & Y < rectangle2_y2) | ...
    (rectangle2_y1 == Y & rectangle2_x1 < X & X < rectangle2_x2) | ...
    (rectangle2_y2 == Y & rectangle2_x1 < X & X < rectangle2_x2));
is_in_border = (is_in_border1 | is_in_border2) & (~is_in_rectangle1_or_2);

figure(1)
plot(X(is_in_rectangle1_or_2),Y(is_in_rectangle1_or_2),'o'); hold on
plot(X(is_in_rectangle1_and_2),Y(is_in_rectangle1_and_2),'o'); hold on
plot(X(is_in_border1),Y(is_in_border1),'o'); hold on
plot(X(is_in_border2),Y(is_in_border2),'o'); hold off

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
figure(2);spy(A_rectangle)
