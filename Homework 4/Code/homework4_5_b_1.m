for k = 1:2
    clearvars -except k; load(strcat("HW4_Problem5b_",int2str(k),".mat"));
    n = length(b); x0 = ones(n,1); tol = 1e-8; k0s = [n,5,10,20];
    clf; m = []; for k0 = k0s
        [~,mm,relres] = gmres_default(A,b,k0,tol,x0);
        m = [m,mm]; plot(log(relres)); hold on;
    end
    xlabel("k"); ylabel("log(\rho_k)"); set(gca,'FontSize',16);
    legend(["No restart",strcat("k_0 = ",string(k0s(2:end)))]);
    title("GMRES (Without Preconditioning)");
    saveas(gcf,strcat("../Figures/homework4_5_",int2str(k),"_default.png"));
    matrix2latex(m',strcat("../Tables/homework4_5_",int2str(k),"_default.tex"),'alignment','r')
    clf; m = []; for k0 = k0s
        [~,mm,relres] = gmres_diag(A,b,k0,tol,x0);
        m = [m,mm]; plot(log(relres)); hold on
    end
    xlabel("k"); ylabel("log(\rho_k)"); set(gca,'FontSize',16);
    legend(["No restart",strcat("k_0 = ",string(k0s(2:end)))]);
    title("GMRES (Diagonal Preconditioning)")
    saveas(gcf,strcat("../Figures/homework4_5_",int2str(k),"_diag.png"));
    matrix2latex(m',strcat("../Tables/homework4_5_",int2str(k),"_diag.tex"),'alignment','r')
    clf; m = []; D = diag(diag(A)); for k0 = k0s
        [~,mm,relres] = gmres_ssor(A,b,k0,tol,D,x0);
        m = [m,mm]; plot(log(relres)); hold on
    end
    xlabel("k"); ylabel("log(\rho_k)"); set(gca,'FontSize',16);
    legend(["No restart",strcat("k_0 = ",string(k0s(2:end)))]);
    title("GMRES (SSOR D = D_0)");
    saveas(gcf,strcat("../Figures/homework4_5_",int2str(k),"_ssordiag.png"));
    matrix2latex(m',strcat("../Tables/homework4_5_",int2str(k),"_ssordiag.tex"),'alignment','r')
    clf; m = []; D = 10*speye(n); for k0 = k0s
        [~,mm,relres] = gmres_ssor(A,b,k0,tol,D,x0);
        m = [m,mm]; plot(log(relres)); hold on
    end
    xlabel("k"); ylabel("log(\rho_k)"); set(gca,'FontSize',16);
    legend(["No restart",strcat("k_0 = ",string(k0s(2:end)))]);
    title("GMRES (SSOR D = 10I)");
    saveas(gcf,strcat("../Figures/homework4_5_",int2str(k),"_ssoriden.png"));
    matrix2latex(m',strcat("../Tables/homework4_5_",int2str(k),"_ssoriden.tex"),'alignment','r')
end