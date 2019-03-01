load('large_ex.mat');
[J,I,V] = comp_l_tri(L);
c = solve_l(J,I,V,b);
[J,I,V] = comp_u_tri(U);
x = solve_u(J,I,V,c);
matrix2latex([x(50000);x(100000);x(150000);x(200000);x(250000)],"../Tables/largeexsolve.tex",'alignment','r','format','%-.15e')