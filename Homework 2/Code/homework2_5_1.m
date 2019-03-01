for k = ["large_ex1","large_ex2"]
    clearvars -except k
    load(strcat(k,".mat"));n=length(b);
    for perm = ["default","colamd","colperm"]
        for scaling = ["default","scaling"]
            if perm == "default"; p0=1:size(A,1);
            else; p0 = eval(strcat(perm,"(A)")); end
            p0i(p0)=1:n;D=speye(n,n);
            if scaling == "default";[L,U,p,q]=lu(A(:,p0),'vector');
            else; [L,U,p,q,D] = lu(A(:,p0),'vector'); end; qi(q)=1:n;
            c = D\b; c = c(p);
            [J,I,V] = comp_l_tri(L); v = solve_l(J,I,V,c);
            [J,I,V] = comp_u_tri(U); x = solve_u(J,I,V,v);
            x = x(qi); x=x(p0i);
            matrix2latex([nnz(L);nnz(U);norm(b-A*x)/norm(b);x(1);x(30000);x(70000);x(140000);x(200002)],strcat("../Tables/",k,"_",scaling,"_",perm,".tex"),'alignment','r','format','%-.15e')
        end
    end
end