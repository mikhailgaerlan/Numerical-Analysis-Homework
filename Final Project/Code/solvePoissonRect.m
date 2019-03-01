function V = solvePoissonRect(F,a,b,c,d)
%Solve two-dimensional Poisson's equation T_m1V + VT_m2 = F
%   input  - F            right-hand side
%            [a,b]        x interval
%            [c,d]        y interval
%   output - V            solution matrix
[m2,m1] = size(F); h1 = (b-a)/(m1+1); h2 = (d-c)/(m2+1);
G = multZ(h2,c,d,multZ(h1,a,b,F.').'); V = zeros(m2,m1);
for k = 1:m1; for j = 1:m2; V(j,k) = G(j,k)/(h1^2*l(h2,j,c,d)+h2^2*l(h1,k,a,b)); end; end
V = multZ(h2,c,d,multZ(h1,a,b,V.').');
end
function lam = l(h,j,a,b); lam = 2*(1-cos(pi*h*j/(b-a))); end
function W = multZ(h,a,b,V)
%Multiply ZV where Z is the eigenvector matrix of T_m1
%   input  - V is an m1 x m2 matrix
%   output - W = ZV
[m1,m2] = size(V);
Vt = [zeros(1,m2);V;zeros(m1+1,m2)];
Wt = fft(Vt); W = -sqrt(2*h/(b-a))*imag(Wt(2:m1+1,:));
end