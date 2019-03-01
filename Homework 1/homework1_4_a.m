load('www0.mat');load('x0.mat');
n = 685230;
y = graphx(n,Q,Jv,x0);
y85 = graphax(n,Q,Jv,0.85,x0);
matrix2latex([y(2);y(222222);y(300000);y(400000)],'www0graphy.tex','alignment','r','format','%-.15e')
matrix2latex([y85(2);y85(222222);y85(300000);y85(400000)],'www0graphy85.tex','alignment','r','format','%-.15e')