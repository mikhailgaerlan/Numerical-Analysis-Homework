load('large_ex.mat');
[J,I,V] = comp_u_tri(U);
matrix2latex([I(50000);I(100000);I(150000);I(200000);I(250000)],"../Tables/largeexu.tex",'alignment','r')