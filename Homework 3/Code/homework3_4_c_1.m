n=10;b=ones(n,1);for p=[1,0.1,0.01]
    t=(1+sqrt(0:(n-1))).^(-p);
    matrix2latex(toeplitzpcg(t,b),strcat("../Tables/pcg_",sprintf('%03d',100*p),".tex"),'alignment','r','format','%-.15e')
    matrix2latex(toeplitzpcgcirc(t,b),strcat("../Tables/pcgcirc_",sprintf('%03d',100*p),".tex"),'alignment','r','format','%-.15e')
end

n=1e6;b=ones(n,1);for p=[1,0.1]
    t=(1+sqrt(0:(n-1))).^(-p);
    x = toeplitzpcg(t,b);
    matrix2latex(real(x([1,100000,500000,700000,1000000])),strcat("../Tables/pcglarge_",sprintf('%03d',100*p),".tex"),'alignment','r','format','%-.15e')
    x = toeplitzpcgcirc(t,b);
    matrix2latex(real(x([1,100000,500000,700000,1000000])),strcat("../Tables/pcgcirclarge_",sprintf('%03d',100*p),".tex"),'alignment','r','format','%-.15e')
end