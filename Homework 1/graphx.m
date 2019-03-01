function y = graphx(n,Q,Jv,x)
%Multiply the matrix representation of the connectivity of a directed graph.
%   input  - n     size of vector
%          - Q     initial connections
%          - Jv    indices of vectors with dj = 0
%          - x     column vector
%   output - A'*x  where A = Q + v*e'/n
A = Q;
y = A'*x + sum(x(Jv'))/n;
end