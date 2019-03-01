load('small_ex.mat');
[J,I,V] = comp_l_tri(L);
c = solve_l(J,I,V,b);
[J,I,V] = comp_u_tri(U);
x = solve_u(J,I,V,c)
matrix2latex(x,"../Tables/smallexsolve.tex",'alignment','r','format','%-.15e')