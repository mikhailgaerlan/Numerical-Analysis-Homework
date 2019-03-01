load('www0.mat');
eps=1e-12;kmax=10000;
n=685230;x10 = ones(n,1);

[lambda,k,releps,~]=powerm(n,Q,Jv,x10,eps,kmax);
matrix2latex([lambda,k,releps]','www0power.tex','alignment','r','format','%-.15e')

[lambda,k,releps,~]=powerm(n,Q,Jv,x0,eps,kmax);
matrix2latex([lambda,k,releps]','www0power0.tex','alignment','r','format','%-.15e')

a = 0.85;
[lambda,k,releps,x]=poweram(n,Q,Jv,a,x10,eps,kmax);
matrix2latex([lambda,k,releps]','www0powera.tex','alignment','r','format','%-.15e')
[~,b]=sort(x,'descend');
matrix2latex(b(1:10)','www0sorta.tex')

[lambda,k,releps,x]=poweram(n,Q,Jv,a,x0,eps,kmax);
matrix2latex([lambda,k,releps]','www0powera0.tex','alignment','r','format','%-.15e')
[~,b]=sort(x,'descend');
matrix2latex(b(1:10)','www0sorta0.tex')