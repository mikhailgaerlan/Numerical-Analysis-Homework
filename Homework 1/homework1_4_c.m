load matrices.mat A;eps=1e-15;kmax=1000;n = 10;
Jv = [3;6;7];Q = A;
for j = Jv'
    Q(j,:) = zeros(1,n);
end

x10 = ones(10,1);
[lambda,k,releps,x]=powerm(n,Q,Jv,x10,eps,kmax);
matrix2latex([lambda,k,releps]','powerm.tex','alignment','r','format','%-.15e')
matrix2latex(x,'powermx.tex','alignment','r','format','%-.15e')
[~,b]=sort(x,'descend');
matrix2latex(b','sortp.tex')

x20 = [1;-1;1;-1;1;-1;1;-1;1;-1];
[lambda,k,releps,x]=powerm(n,Q,Jv,x20,eps,kmax);
matrix2latex([lambda,k,releps]','powerm2.tex','alignment','r','format','%-.15e')
matrix2latex(x,'powerm2x.tex','alignment','r','format','%-.15e')
[~,b]=sort(x,'descend');
matrix2latex(b','sortp2.tex')