m = 10;[X,Y] = meshgrid(1:m);
h = 1/(m+1);l(h,


function lam = l(h,j)
    lam = 2.*(1-cos(pi.*h.*j));
end