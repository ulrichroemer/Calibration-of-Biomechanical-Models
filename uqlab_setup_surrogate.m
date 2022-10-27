function [] = uqlab_setup_surrogate()

    PCE_opts.ExpDesign.X = theta_ed;
    PCE_opts.ExpDesign.Y = xi;
    PCE_opts.Type = 'Metamodel';    
    PCE_opts.MetaType = 'PCE';            
    PCE_opts.Method = 'OLS';

    for i = 1:M
        IOpts.Marginals(i).Type = 'Uniform';
        IOpts.Marginals(i).Parameters = [-5,1];    
    end
    myInput = uq_createInput(IOpts);
    PCE_opts.Degree = P;       
    PCE_metamodel = uq_createModel(PCE_opts);   


end