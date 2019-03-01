clear variables; load('pcg_small.mat'); global J I V; n = size(A,1);
[x,~,~,iter] = pcg(A,b,1e-9,1000,[],[],ones(n,1));
matrix2latex(x,"../Tables/snocond.tex",'alignment','r','format','%-.15e');
fileID = fopen('../Tables/snoconditer.tex','w'); fprintf(fileID,'%d',iter);fclose(fileID);
fileID = fopen('../Tables/snocondnorm.tex','w'); fprintf(fileID,'%.15e',norm(A*x-b)/norm(b));fclose(fileID);

[J,I]=get_lower(A); V = iCholesky(A,J,I); col = 3; m = length(J);
rows = ceil(m/col); gap = rows*col-m;
[x,~,~,iter] = pcg(A,b,1e-9,10000,@solve_lowerr,@solve_lowertt,ones(n,1));
matrix2latex(x,"../Tables/scond.tex",'alignment','r','format','%-.15e');
fileID = fopen('../Tables/sconditer.tex','w'); fprintf(fileID,'%d',iter);fclose(fileID);
fileID = fopen('../Tables/scondnorm.tex','w'); fprintf(fileID,'%.15e',norm(A*x-b)/norm(b));fclose(fileID);
for j = (m+1):(m+gap); J(j) = NaN; V(j) = NaN; end
matrix2latex(reshape(J,[rows,col]),"../Tables/scondJ.tex",'alignment','r','format','%-d');
matrix2latex(reshape(I,[13,2]),"../Tables/scondI.tex",'alignment','r','format','%-d');
matrix2latex(reshape(V,[rows,col]),"../Tables/scondV.tex",'alignment','r','format','%-.15e');

load('pcg_large.mat'); n = size(A,1);
[x,~,~,iter] = pcg(A,b,1e-9,10000,[],[],ones(n,1));
matrix2latex(x([1,10000,100000,20000,262144]),"../Tables/lnocond.tex",'alignment','r','format','%-.15e');
fileID = fopen('../Tables/lnoconditer.tex','w'); fprintf(fileID,'%d',iter);fclose(fileID);

[J,I]=get_lower(A); V = iCholesky(A,J,I);
[x,~,~,iter] = pcg(A,b,1e-9,10000,@solve_lowerr,@solve_lowertt,ones(n,1));
matrix2latex(x([1,10000,100000,20000,262144]),"../Tables/lcond.tex",'alignment','r','format','%-.15e');
fileID = fopen('../Tables/lconditer.tex','w'); fprintf(fileID,'%d',iter);fclose(fileID);
matrix2latex(J([1,10,100,1000,end]),"../Tables/lcondJ.tex",'alignment','r','format','%-d');
matrix2latex(I([1,10,100,1000,end]),"../Tables/lcondI.tex",'alignment','r','format','%-d');
matrix2latex(V([1,10,100,1000,end]),"../Tables/lcondVL.tex",'alignment','r','format','%-.15e');

function y=solve_lowerr(xx); global J I V; y=solve_lower(xx,J,I,V); end
function y=solve_lowertt(xx); global J I V; y=solve_lowert(xx,J,I,V); end