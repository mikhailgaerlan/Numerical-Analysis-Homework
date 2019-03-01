load('small_ex.mat');A=L*U;
p0=colamd(A);
[L,U,p,q,D] = lu(A(:,p0),'vector');
b=D\b;b=b(p);
[J,I,V] = comp_l_tri(L);
a = solve_l(J,I,V,b);
[J,I,V] = comp_u_tri(U);
a = solve_u(J,I,V,a);
for i = 1:length(b); y=L*U*a; fprintf("%25.15e %25.15e\n", y(i),b(i)); end
