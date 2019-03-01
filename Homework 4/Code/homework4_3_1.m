clearvars; clf; load('HW4_Problem3.mat');
n = length(b); k = [n,2,5,10,20,50,100]; m = []; x0 = zeros(n,1);
for k0 = k
    [~,n,relres] = gmres_default(A,b,k0,1e-9,x0);
    m = [m,n]; plot(log(relres)); hold on;
end
xlabel("k"); ylabel("log(\rho_k)"); set(gca,'FontSize',16);
title("GMRES"); legend(["No restart",strcat("k_0 = ",string(k(2:end)))]);
saveas(gcf,"../Figures/homework4_3.png");
matrix2latex(m',"../Tables/homework4_3.tex",'alignment','r')