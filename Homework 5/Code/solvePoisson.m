function V = solvePoisson(m,F)
%Solve two-dimensional Poisson's equation T_mV + VT_m = F
%   input  - m            matrix size
%            F            right-hand side
%   output - V            solution matrix
h = 1/(m+1); G = multZ(h,multZ(h,F.').'); V = zeros(m,m);
for k = 1:m; for j = 1:m; V(j,k) = G(j,k)/(l(h,j)+l(h,k)); end; end
V = multZ(h,multZ(h,V.').');
end
function lam = l(h,j); lam = 2*(1-cos(pi*h*j)); end
function W = multZ(h,V)
%Multiply ZV where Z is the eigenvector matrix of T_m
%   input  - V is an n x n matrix
%   output - W = ZV
m = size(V,1);
Vt = vertcat(zeros(1,m),V,zeros(m+1,m));
Wt = fft(Vt);
W = -sqrt(2*h)*imag(Wt(2:m+1,:));
end