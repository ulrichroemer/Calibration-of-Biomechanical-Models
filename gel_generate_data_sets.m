%% Generate LHS sample of input parmaters and store

N_MC = 1000;
n = 3;    % number of parameters

% lin-scale
% p_min = 1e-5; 
% p_max = 1e1; 
% p = @(s) p_min + (p_max-p_min)*(1+s)/2;

% log-scale
p = @(s) 10.^(s);

S = -5 + 6*lhsdesign(N_MC,n);
P = p(S);

%% generate id
temp =  java.util.UUID.randomUUID;
myuuid = temp.toString

fid = strcat(['data/gel_exp_design_input_material_id_', char(myuuid), '.csv']);
csvwrite(fid,P)
