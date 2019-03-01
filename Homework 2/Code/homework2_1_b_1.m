m = 37;nnzs = zeros(5,1);A = make_3d_laplacian(m);

orderings = ["default","symamd","colamd","symrcm","colperm"];n=1;
for order = orderings
    if n == 1; p = 1:m^2;
    else; p = eval(strcat(order,"(A)")); end
    L = chol(A(p,p),'lower');
    nnzs(n) = nnz(L);spy(L);
    if n == 1; title("No reordering",'FontName','Times');
    else; title(order,'FontName','Courier'); end; set(gca,'FontSize',20);
    saveas(gcf,strcat("../Figures/",order,"b"),'png');n=n+1;
end; matrix2latex(nnzs,"../Tables/nnzsb.tex",'alignment','r')