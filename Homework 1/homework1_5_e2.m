load('t_and_x.mat')
y = toeplitz(t,x);
matrix2latex([y(1),y(100000),y(200000),y(300000),y(400000),y(500000),sum(y),norm(y,2)]','t_and_x.tex','alignment','r','format','%-.15e')