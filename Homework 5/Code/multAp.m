function v = multAp(x,gamma,m)
%Multiply a vector v = (A_0 + gamma A_1) x
%   input  - x      input vector
%          - gamma  constant
%          - m      size of vector x
%   output - v      output vector
h = 1/(m+1); v = x + gamma*solveM1(multA1(x));
    function v = solveM1(x)
        H = reshape(x,[m,m]);
        V = solvePoisson(m,H);
        v = reshape(V,[m*m,1]);
    end
    function v = multA1(b)
        v = vertcat(b(m+1:end),zeros(m,1));
        v(m+1:end) = v(m+1:end)-b(1:end-m);
        v = v*h/2;
    end
end