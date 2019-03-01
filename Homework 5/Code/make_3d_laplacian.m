function A = make_3d_laplacian(m)

n = m^3;

m2 = m^2;

B = zeros(n,4);

B(:,1) = 3 * ones(n,1);

Xtmp = [2 : m]' * ones(1,m2) + m * ones(m-1,1) * [0 : m2-1];
I = reshape(Xtmp, n-m2, 1);
B(I,2) = -ones(n-m2,1);

Xtmp = [m+1 : m2]' * ones(1,m) + m2 * ones(m2-m,1) * [0 : m-1];
I = reshape(Xtmp, n-m2, 1);
B(I,3) = -ones(n-m2,1);

B(:,4) = -ones(n,1);

A = spdiags(B,[0 1 m m2],n,n);

A = A + A';

