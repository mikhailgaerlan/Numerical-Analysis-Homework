format long e;
t = [6,-5,1,2,-1,4,-3];
T = [2,-1,4,-3;...
    1,2,-1,4;...
    -5,1,2,-1;...
    6,-5,1,2];
x = [-1,2,1,4]';
matrix2latex(T*x,'Tx.tex','alignment','r','format','%-.15e')
matrix2latex(toeplitz(t,x),'toeplitzx.tex','alignment','r','format','%-.15e')