clearvars; m = 100; h = 1/(m+1); [Y,X]=meshgrid(h:h:1-h);
tol = 1e-10; ms = []; relress = [];
x = X(1,:);y = Y(:,1); F = f(X,Y);
b0 = ones(1,m); b1 = ones(1,m); c0 = ones(m,1); c1 = ones(m,1);
for gamma = [1,10,50,100,1000]
    G = h^2*F; G(1,:) = G(1,:)+b0; G(end,:) = G(end,:)+b1;
    G(:,1) = G(:,1)+c0+gamma*c0*h/2; G(:,end) = G(:,end)+c1-gamma*c1*h/2;
    [v,iter,relres] = poissonSolve(m,reshape(G,[m*m,1]),gamma,tol);
    ms = [ms,iter]; relress = [relress,relres]; clf; colorbar;
    surf(X,Y,reshape(v,[m,m]),'LineStyle','none'); hold on; view(2);
    title(strcat("\gamma = ",int2str(gamma))); set(gca,'FontSize',20);
    saveas(gcf,strcat("../Figures/homework5_1_",int2str(gamma),".png"));
end
matrix2latex(ms','../Tables/homework5_1_iters.tex','alignment','r')
matrix2latex(relress','../Tables/homework5_1_relres.tex','alignment','r','format','%-.15e')
function fun = f(x,y); fun = x.^3.*y.^2.*exp(2-x-y); end