format long e;
load matrices.mat A;
[V,D]=eig(A');l = diag(D);
[absl,c]=sort(abs(l));
x = abs(V(:,c(end):c(end)));
[~,b]=sort(x,'descend');
matrix2latex(x,'x.tex','format','%-.15e')
matrix2latex(b','sort.tex')