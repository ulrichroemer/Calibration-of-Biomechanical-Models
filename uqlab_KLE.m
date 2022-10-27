function [lambda,Vs,Es]= uqlab_KLE(G_ed,T)

    [V,E] = eig(cov(G_ed));
    [e,ind] = sort(diag(E),'descend');
    Vs = V(:,ind);
    Es = E(ind,ind);
    lambda = diag(Es);
    lambda = lambda(1:T);
    Vs = Vs(:,1:T);
    Es = Es(1:T,1:T);
    
    
end

