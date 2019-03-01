load("large_ex2.mat");n=length(b);
for perm = ["default","colamd","colperm"]
    for scaling = ["default","scaling"]
        
        if perm == "default"; p0=1:size(A,1);
        else; p0 = eval(strcat(perm,"(A)")); end
        p0i(p0)=1:n;D=speye(n);
        if scaling == "default";[L,U,p,q]=lu(A(:,p0),'vector');
        else; [L,U,p,q,D] = lu(A(:,p0),'vector'); end; qi(q)=1:n;
        
        c = D\b; c = c(p);
        
        [J,I,V] = comp_l_tri(L); v = solve_l(J,I,V,c);
        [J,I,V] = comp_u_tri(U); x = solve_u(J,I,V,v);
        
        x = x(qi); x=x(p0i);
        %for i = 1:length(b); fprintf("%25.15e %25.15e\n", x(i),y(i)); end
        norm(b-A*x)/norm(b)
    end
end
