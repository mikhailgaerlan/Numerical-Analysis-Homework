m0=37;i=0;
while 1    
    fileID = fopen('../Tables/ib.tex','w');
    fprintf(fileID,int2str(i));fclose(fileID);
    m=2^i*m0;A = make_3d_laplacian(m);
    p = symamd(A);
    L = chol(A(p,p),'lower');
    i = i + 1;
end