clearvars; clc; m = 8; nmax = 100; tol = 1e-6;
tic; params = [1/2,0,0,1,1/2;1/2,0,2,7/2,2;1/4,1,0,1,1/2;1/4,1,2,7/2,2;...
    1/2,2,0,1,1/2;1/2,2,2,7/2,2]; n = size(params,1);
relerrs1 = []; relerrs2 = []; nsolves1 = []; nsolves2 = [];
for i = 1:n
    as=params(i,1);bs=params(i,2);alpha=params(i,3);beta=params(i,4);
    gamma=params(i,5);
    [A,~,~,I1,I2,X,Y,b] = laplace_matching(m,as,bs,...
        @(x,y) f(x,y,alpha,beta,gamma),@(x,y) v(x,y,alpha,beta,gamma));
    V = v(X,Y,alpha,beta,gamma);
    [Va1,iters1] = mult_Schwarz(m,nmax,as,bs,tol,A,b,I1,I2);
    [Va2,iters2] = CG_add_Schwarz(m,nmax,as,bs,tol,A,b,I1,I2);
    err1=max(abs(Va1-V))/max(abs(V));err2=max(abs(Va2-V))/max(abs(V));
    relerrs1 = [relerrs1;err1]; relerrs2 = [relerrs2;err2];
    nsolves1 =[nsolves1;iters1]; nsolves2 =[nsolves2;iters2];
    scatter3(X,Y,V,'o','fill','CData',V);set(gca,'FontSize',16);
    title(strcat("Actual Case ",int2str(i)));
    saveas(gcf,strcat("../Figures/final_3_actual_",int2str(i),".png"));
    scatter3(X,Y,Va1,'o','fill','CData',Va1);set(gca,'FontSize',16);
    title(strcat("Multiplicative Schwarz Case ",int2str(i)));
    saveas(gcf,strcat("../Figures/final_3_mult_",int2str(i),".png"));
    scatter3(X,Y,Va2,'o','fill','CData',Va2);set(gca,'FontSize',16);
    title(strcat("CG with Additive Schwarz Case ",int2str(i)));
    saveas(gcf,strcat("../Figures/final_3_pcg_",int2str(i),".png"));
end
matrix2latex(relerrs1,'../Tables/final_3_relerr_mult.tex','alignment','r','format','%-.15e');
matrix2latex(relerrs2,'../Tables/final_3_relerr_pcg.tex','alignment','r','format','%-.15e');
matrix2latex(nsolves1,'../Tables/final_3_iters_mult.tex','alignment','r','format','%-4d');
matrix2latex(nsolves2,'../Tables/final_3_iters_pcg.tex','alignment','r','format','%-4d');toc;
function true = v(x,y,a,b,g); true=y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y); end
function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end