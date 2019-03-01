clearvars; load("HW5_P5b.mat"); m = 69;
A = make_3d_laplacian(m);
for kmax = [250,500,1000,2000]
    clearvars -except kmax A m r;
    Tk = hermlanc(A,r,kmax,1e-14); eigs = eig(Tk);
    for i=1:m;for j=i:m; for l=j:m;ls(i,j,l)=le(i,j,l,m);end;end;end
    ls = ls(ls~=0);
    matrix2latex(mink(eigs,10),strcat("../Tables/homework5_5_b_min_approx_",num2str(kmax),".tex"),'alignment','r','format','%-.15e')
    matrix2latex(mink(ls,10),strcat("../Tables/homework5_5_b_min_exact_",num2str(kmax),".tex"),'alignment','r','format','%-.15e')
    matrix2latex(maxk(eigs,10),strcat("../Tables/homework5_5_b_max_approx_",num2str(kmax),".tex"),'alignment','r','format','%-.15e')
    matrix2latex(maxk(ls,10),strcat("../Tables/homework5_5_b_max_exact_",num2str(kmax),".tex"),'alignment','r','format','%-.15e')
end
function lam = le(i,j,l,m)
lam = 2.*(3-cos((i.*pi)./(m+1))-cos((j.*pi)./(m+1))-cos((l.*pi)./(m+1)));
end