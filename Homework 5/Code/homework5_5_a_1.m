clearvars; load("HW5_P5a.mat"); m = 2;
A = make_3d_laplacian(m); kmax = 2000; Tk = hermlanc(A,r,kmax,1e-10);
matrix2latex(eig(Tk),"../Tables/homework5_5_a_approx.tex",'alignment','r','format','%-.15e')
for i=1:m;for j=i:m; for l=j:m;lambdas(i,j,l)=lambda(i,j,l,m);end;end;end
matrix2latex(lambdas(lambdas~=0),"../Tables/homework5_5_a_exact.tex",'alignment','r','format','%-.15e')
function lam = lambda(i,j,l,m)
lam = 2.*(3-cos((i.*pi)./(m+1))-cos((j.*pi)./(m+1))-cos((l.*pi)./(m+1)));
end