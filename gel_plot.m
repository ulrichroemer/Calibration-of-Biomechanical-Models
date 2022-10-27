%% Plotting gel results
clearvars

addpath(genpath(pwd))
close all
rng('default')
verbose = 1;
export = 1;

%n = 1;    % number of parameters
%[x_sim,F_sim_1] = gel_solver(0.01,n);    % get correct loading values

%n = 2;    % number of parameters
%[~,F_sim_2] = gel_solver([0.01,0.001],n);    % get correct loading values

n = 3;    % number of parameters
[x_sim,F_sim_3] = gel_solver([0.01,0.01, 0.01],n);    % get correct loading values

% load measurement data
[x_meas,f_meas] = get_data();
b = interp1(x_meas,f_meas,x_sim);
b = b(2:end)';
b = log(b);

% plot simulated vs measured data
%figure
%plot(x_sim(2:end),b,'r--',x_sim,F_sim_1,'k+-',x_sim,F_sim_2,'k*-',x_sim,F_sim_3,'ko-')

%% Setup
%id = 'a7ea9e74-48be-47a6-847e-cf99a66bea29';    % 100 (n=1); log-scale 
%id = '06c29bfc-9b12-4a17-9a7e-8f38cc24f8a7';    % 100 (n=2); log-scale 
%id = 'f4e35095-7b8b-4a55-bf4c-845939eade29';    % 100 (n=3); log-scale  
id = '0c3f341d-2682-493c-9674-43cc5c0c6b4e';    % 1000 (n=3); log-scale  

id_in = ['gel_exp_design_input_material_id_',id,'.csv'];
id_out = ['gel_simulated_data_id_',id,'.csv'];

ind = 900;    % index, where experimental design is split

F_ed = csvread(id_out);
theta_ed = csvread(id_in);

F_cv = F_ed(2:end,(ind+2):end)';    % first column contains displacements
F_ed = F_ed(2:end,2:(ind+1))';

theta_cv = theta_ed((ind+1):end,:);
theta_ed = theta_ed(1:ind,:);

theta_cv = rescale(theta_cv);
theta_ed = rescale(theta_ed);

% apply log-transformation
G_cv = log(F_cv(:,1:end));
G_ed = log(F_ed(:,1:end));

K = size(G_ed,1);    % number of samples
K_cv = size(G_cv,1);    % number of samples
M = size(theta_cv,2);    % number of parameters

rng('default')
f = figure;
plot(x_sim(2:end),G_cv(randi(100,30,1),:)')
%title('Force responses in the experimental design')
%xlabel('loading step i')
%ylabel('G_i')
if export 
    set(gcf,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
    %set(h, 'Position', [0.2, 0.8, 0.05, .1])
    print(f,'Exp_design.png','-dpng','-r400') 
end

%% Plot KLE 
T = 20;
[lambda,Vs] = KLE(G_ed,T);

% plot eigenvalues
close all
f = figure;
semilogy(1:20,lambda(1:20)/lambda(1),'*')        
if export
    set(f,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
    print(f,'Eigenvalues.png','-dpng','-r400') 
end

% plot eigenfunctions
g = figure;
ii = 1:length(Vs(:,1));
plot(ii,sqrt(lambda(1))*Vs(:,1),ii,sqrt(lambda(2))*Vs(:,2),ii,sqrt(lambda(3))*Vs(:,3),ii,sqrt(lambda(4))*Vs(:,4))        
h = legend({'$\sqrt{\lambda_1}\mathbf{v}_1$','$\sqrt{\lambda_2}\mathbf{v}_2$','$\sqrt{\lambda_3}\mathbf{v}_3$','$\sqrt{\lambda_4}\mathbf{v}_4$'},'interpreter','latex');
if export
    set(g,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
    set(h, 'Position', [0.3, 0.65, 0.05, .1])
    print(g,'Eigenfunctions.png','-dpng','-r400') 
end

%% Plot surrogate against ED
T = 3;
P = 5;
mdl = setup_surrogate(M,T,P,G_ed,theta_ed);
[lambda,Vs,Es] = KLE(G_ed,T);
G_surr = @(t) (mean(G_ed)' + Vs(:,1:T)*Es(1:T,1:T).^(1/2)*uq_evalModel(mdl,t)')';

f = figure 
plot(x_sim(2:end),b,'r+')
hold on
C = {'k','b','m','g'};
for i = 1:4
    delta = 1;
    plot(x_sim(2:end),G_cv(i+delta,:),'k-',x_sim(2:end),G_surr(theta_cv(i+delta,:)),'k--')
end
h = legend('data','model $G_i^{(k)}$', 'surrogate $\tilde{G}_i^{(k)}$', 'Interpreter', 'latex');
if export
    set(f,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
    set(h, 'Position', [0.7, 0.2, 0.05, .1])
    print(f,'Surrogate_vs_Design.png','-dpng','-r400') 
end

%% Plot surrogate accuracy 

T_max = 5;
P_max = 5;
error = zeros(T_max,P_max);

for i = 1:T_max
    for j = 1:P_max
        mdl = setup_surrogate(M,i,j,G_ed,theta_ed);
        [lambda,Vs,Es] = KLE(G_ed,i);
        G_surr = @(t) (mean(G_ed)' + Vs(:,1:i)*Es(1:i,1:i).^(1/2)*uq_evalModel(mdl,t)')';
        error(i,j) = norm(G_surr(theta_cv)-G_cv);
    end
end

%%
close all
f = figure;
semilogy(1:T_max,error(:,end),'*-')        
if export
    set(f,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
    print(f,'Error_KLE.png','-dpng','-r400') 
end

close all
f = figure;
semilogy(1:P_max,error(end,:),'*-')        
if export
    set(f,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
    print(f,'Error_PCE.png','-dpng','-r400') 
end


%% aux functions

function [theta_rs] = rescale(theta)

    % lin-scale
%     pmin = 0.001; 
%     pmax = 0.5;
%     theta_rs = 2*((theta - pmin)/(pmax - pmin) - 1) + 1;
    
    % log-scale
    theta_rs = log10(theta);
    
end

function [theta] = upscale(theta_rs)

    % lin-scale
%     pmin = 0.001; 
%     pmax = 0.5;
%     theta(:,1) = (theta_rs + 1)*(pmax - pmin)/2 + pmin;
    
    % log-scale
    theta = 10.^(theta_rs);    
end

function [x_meas,f_meas,sigma_meas] = get_data()

    % get the exp data
    fpath = cd;
    exppath = 'ISM_library\Liu\exp';
    cd(exppath)
    MD_all = xlsread('exp.xlsx');
    MD1 = [MD_all(:,1),MD_all(:,2)]; % strain and force
    cd(fpath)
    x_meas = MD_all(:,1);
    f_meas = mean(MD_all(:,2:end),2);
    sigma_meas = std(MD_all(:,2:end)')';
    %f_meas = MD_all(:,2);
    
    x_meas = x_meas(1:6:end);
    f_meas = f_meas(1:6:end);
    sigma_meas = sigma_meas(1:6:end);
    
end

function [PCE_metamodel] = setup_surrogate(M,T,P,G_ed,theta_ed)

    
    %explained_var = sum(lambda(1:T))/sum(lambda)

    [lambda,Vs] = KLE(G_ed,T);
    xi = (diag(1./sqrt(lambda))*Vs'*(G_ed - mean(G_ed))')';

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
    %uq_display(PCE_metamodel,1)

    
end

function [lambda,Vs,Es]= KLE(G_ed,T)

    [V,E] = eig(cov(G_ed));
    [e,ind] = sort(diag(E),'descend');
    Vs = V(:,ind);
    Es = E(ind,ind);
    lambda = diag(Es);
    lambda = lambda(1:T);
    Vs = Vs(:,1:T);
    
end


