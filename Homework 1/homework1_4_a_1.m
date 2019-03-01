load 'matrices.mat' A
n = 10;
Jv = [3;6;7];
Q = A;
for j = Jv'
    Q(j,:) = zeros(1,n);
end
matrix2latex(graphx(n,Q,Jv,ee),'graphy.tex','format','%-.15e')
matrix2latex(graphax(n,Q,Jv,0.5,ee),'graphy05.tex','format','%-.15e')
matrix2latex(graphax(n,Q,Jv,0.85,ee),'graphy085.tex','format','%-.15e')