clearvars; m = 100; h = 1/(m+1); kmax = 300;
load("HW5_P4.mat"); mins=[]; maxes=[]; minlambdas=[]; maxlambdas=[];
for gamma = [1,10,50,100,1000]
    [Hkt,Vkk] = arnoldi(@(z)multAp(z,gamma,m),r,kmax,0);
    Hk = Hkt(1:end-1,:); [Z,Lt] = eig(Hk); lambdas = diag(Lt);
    rhos = Hkt(end,end)*abs(Z(end,:));
    clf; plot(lambdas,'.');
    title(strcat("\gamma = ",int2str(gamma))); set(gca,'FontSize',20);
    saveas(gcf,strcat("../Figures/homework5_4_",int2str(gamma),".png"));
    rhomin = min(rhos); rhomax = max(rhos);
    mins = [mins;rhomin]; maxes = [maxes;rhomax];
    i = find(rhos==rhomin); j = find(rhos==rhomax);
    matrix2latex(lambdas(i),strcat("../Tables/homework5_4_mins_",num2str(gamma),".tex"),'alignment','r','format','%-.15e')
    matrix2latex(lambdas(j),strcat("../Tables/homework5_4_maxes_",num2str(gamma),".tex"),'alignment','r','format','%-.15e')
end
matrix2latex(mins,"../Tables/homework5_4_mins.tex",'alignment','r','format','%-.15e')
matrix2latex(maxes,"../Tables/homework5_4_maxes.tex",'alignment','r','format','%-.15e')