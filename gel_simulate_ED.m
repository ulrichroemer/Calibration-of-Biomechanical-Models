% Simulate experimental design

%% Setup
addpath(genpath(pwd))
clc; close all; clear all

% Call solver 
n = 3;    % number of parameters
pset = [0.003,0.003,0.003];    % unit [MPa]
[e_fixer,f_fixer] = gel_solver(pset,n);
n_load = length(f_fixer);

%% LHS sample 
id = '0c3f341d-2682-493c-9674-43cc5c0c6b4e';
id_in = ['gel_exp_design_input_material_id_',id,'.csv'];
id_out = ['gel_simulated_data_id_',id,'.csv'];
p_ED = csvread(['data\',id_in]);
n_ED = length(p_ED);
f = zeros(n_load,n_ED);

for i = 1:n_ED
    i
    [e_loc,f_fixer] = gel_solver(p_ED(i,:),n);
    if length(f_fixer) ~= length(e_fixer)
        f_fixer = interp1(e_loc,f_fixer,e_fixer);
    end
    f(:,i) = f_fixer;
end

csvwrite(['data\',id_out],[e_fixer,f])