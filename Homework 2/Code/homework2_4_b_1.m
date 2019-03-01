load('small_ex.mat');
[J,I,V] = comp_u_tri(U);
matrix2latex(J,"../Tables/smallexuj.tex",'alignment','r')
matrix2latex(I,"../Tables/smallexui.tex",'alignment','r')
matrix2latex(V,"../Tables/smallexuv.tex",'alignment','r','format','%-.15e')