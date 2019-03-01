load('small_ex.mat');
[J,I,V] = comp_l_tri(L);
matrix2latex(J,"../Tables/smallexj.tex",'alignment','r')
matrix2latex(I,"../Tables/smallexi.tex",'alignment','r')
matrix2latex(V,"../Tables/smallexv.tex",'alignment','r','format','%-.15e')