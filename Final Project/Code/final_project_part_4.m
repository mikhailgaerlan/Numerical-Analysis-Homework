clearvars; clc; m = 7; nmax = 30; tol = 1e-6;
tic; params = [2/5,1/5,0,1,1/2;2/5,1/5,2,7/2,2;...
    1/5,6/5,0,1,1/2;1/5,6/5,2,7/2,2;...
    2/5,9/5,0,1,1/2;2/5,9/5,2,7/2,2]; n = size(params,1);
relerrs1 = []; relerrs2 = [0]; nsolves1 = []; nsolves2 = [0];
for i = 1:n
    as=params(i,1);bs=params(i,2);alpha=params(i,3);beta=params(i,4);
    gamma=params(i,5);
    [vapp,X,Y,iters1] = alt_Schwarz(m,nmax,as,bs,tol,@(x,y)f(x,y,alpha,beta,gamma),@(x,y)v(x,y,alpha,beta,gamma));
    vact = v(X,Y,alpha,beta,gamma);
    err1=max(abs(vact-vapp))/max(abs(vact));
    relerrs1 = [relerrs1;err1]; nsolves1 =[nsolves1;iters1];
    scatter3(X,Y,vact,'o','fill','CData',vact);set(gca,'FontSize',16);
    title(strcat("Actual Case ",int2str(i)));
    saveas(gcf,strcat("../Figures/final_4_actual_",int2str(i),".png"));
    scatter3(X,Y,vapp,'o','fill','CData',vapp);set(gca,'FontSize',16);
    title(strcat("Discretized Alternating Schwarz Case ",int2str(i)));
    saveas(gcf,strcat("../Figures/final_4_add_",int2str(i),".png"));
end
matrix2latex(relerrs1,'../Tables/final_4_relerr_add.tex','alignment','r','format','%-.15e');
matrix2latex(relerrs2,'../Tables/final_4_relerr_gmres.tex','alignment','r','format','%-.15e');
matrix2latex(nsolves1,'../Tables/final_4_iters_add.tex','alignment','r','format','%-4d');
matrix2latex(nsolves2,'../Tables/final_4_iters_gmres.tex','alignment','r','format','%-4d'); toc;
function true = v(x,y,a,b,g); true=y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y); end
function fun = f(x,y,a,b,g)
fun = b.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y)...
    -a.*(a-1).*y.^(a-2).*sin(b.*pi.*x).*cos(g.*pi.*y)...
    +2.*a.*g.*pi.*y.^(a-1).*sin(b.*pi.*x).*sin(g.*pi.*y)...
    +g.^2.*pi.^2.*y.^a.*sin(b.*pi.*x).*cos(g.*pi.*y);
end