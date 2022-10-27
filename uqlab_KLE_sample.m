function [xi_design]= uqlab_KLE_sample(eigenvalues,eigenvectorMatrix,Y_design)
    
    xi_design = (diag(1./sqrt(eigenvalues))*eigenvectorMatrix'*(Y_design - mean(Y_design))')';
    
end

