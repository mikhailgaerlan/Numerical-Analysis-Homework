function [vapp,k] = mult_Schwarz(m,nmax,as,bs,tol,A,b,I1,I2)
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

vapp = zeros(length(b),1);
error = b - A*vapp;
for n = 1:nmax
    vapp = vapp + multB1(b-A*vapp);
    vapp = vapp + multB2(b-A*vapp);
    error = error-Minv(A*error); k = n; if norm(error) < tol; break; end
end

    function b = multB1(v)
        m1 = rectangle1_m1-2; m2 = rectangle1_m2-2;
        z = reshape(I1*v,[m2,m1]);
        g = solvePoisson(z,rectangle1_x1,rectangle1_x2,...
            rectangle1_y1,rectangle1_y2);
        b = I1'*reshape(g,[m1*m2,1]);
    end
    function b = multB2(v)
        m1 = rectangle2_m1-2; m2 = rectangle2_m2-2;
        z = reshape(I2*v,[m2,m1]);
        g = solvePoisson(z,rectangle2_x1,rectangle2_x2,...
            rectangle2_y1,rectangle2_y2);
        b = I2'*reshape(g,[m1*m2,1]);
    end
    function b = Minv(v); b = multB1(v)+multB2(v)-multB2(A*multB1(v)); end
end