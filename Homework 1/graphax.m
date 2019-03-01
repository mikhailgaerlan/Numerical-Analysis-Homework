function y = graphax(n,Q,Jv,a,x)
%Multiply the matrix representation of the connectivity of a directed graph.
%   input  - n     size of vector
%          - Q     initial connections
%          - Jv    indices of vectors with dj = 0
%          - a     parameter between 0 and 1
%          - x     column vector
%   output - A'*x  where A = Q + v*e'/n
A = Q;
y = a*A'*x + a*sum(x(Jv'))/n + (1-a)*sum(x)/n;
end