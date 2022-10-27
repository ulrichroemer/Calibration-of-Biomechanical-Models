function [e_fixer,f_fixer] = gel_solver(pset,n)

if(length(pset)~=n)
       disp('Parameter number does not match n')
       return
end

%% define the initial parameter set

cd('C:\Users\URoemer\bayesian-calibration-of-oocytes');

simpath = 'ISM_library\Liu\sim\';

%% write the parameter and run simulation with the rewritten input file

cd(simpath)
qq = 1;
switch n 
    
    case 1 
        origfile = 'WPG_Zug_M1.inp';
        inputfile = ['WPG_Zug_M1' num2str(qq) '.inp'];   
        outputfile = ['cWPG_Zug_M1' num2str(qq) '.inp'];
    case 2
        origfile= ['WPG_Zug_M2.inp'];
        inputfile= ['WPG_Zug_M2' num2str(qq) '.inp'];   
        outputfile= ['cWPG_Zug_' num2str(qq) '.inp'];
    case 3
        origfile= ['WPG_Zug_M3.inp'];
        inputfile= ['WPG_Zug_M3' num2str(qq) '.inp'];   
        outputfile= ['cWPG_Zug_' num2str(qq) '.inp'];
end


data = importdata(origfile,'\n'); 
%find the parameters of ogden model
row_idx = find(~cellfun('isempty',strfind(data,'Hyperelastic')));
paraline = row_idx+1;
fid = fopen(origfile);
txt = textscan(fid,'%s','delimiter','\n');
fclose(fid);
rows = txt{1,1};
% hyperelastic
M = pset;                   %parameter matrix

switch n 
    case 1 
        formatsp='%d,   0';    % parameter format
        %   0.001, 0.002,    0.,    0.
        M0 = {sprintf(formatsp,M(1))};   % transform
    case 2 
        formatsp='%d,  %d, 0.,0.';    % parameter format
        M0= {sprintf(formatsp,M(1),M(2))};   % transform
    case 3
        formatsp='%d,  %d, %d, 0., 0., 0.';    % parameter format
        M0= {sprintf(formatsp,M(1),M(2),M(3))};   % transform
end

A = rows(1:paraline-1);
B = rows(paraline+1:end);
rows = [A;M0;B];
fid = fopen(outputfile,'wt');
fprintf(fid, '%s\n', rows{:});
fclose(fid);
delete(inputfile)

copyfile(outputfile,inputfile);


%qq %show the set number of the generation
switch n 
    case 1 
        vorhanden=exist('WPG_sim_M1.lck');         %%%%%%//----------------\\%
        if vorhanden==2                         % <===   check the lock   %
            delete('WPG_sim_M1.lck');               %%%%%%\\----------------//%
        end
    case 2
        vorhanden=exist('WPG_sim_M2.lck');         %%%%%%//----------------\\%
         if vorhanden==2                         % <===   check the lock   %
            delete('WPG_sim_M2.lck');               %%%%%%\\----------------//%
         end     
    case 3 
        vorhanden=exist('WPG_sim_M3.lck');         %%%%%%//----------------\\%
         if vorhanden==2                         % <===   check the lock   %
            delete('WPG_sim_M3.lck');               %%%%%%\\----------------//%
         end     
end


%% catch the simulated results

switch n 
    case 1 
        cmd_str = ['abaqus job=WPG_sim_M1 input=WPG_Zug_M1', num2str(qq), ' interactive'];
    case 2 
        cmd_str = ['abaqus job=WPG_sim_M2 input=WPG_Zug_M2', num2str(qq), ' interactive'];
    case 3
        cmd_str = ['abaqus job=WPG_sim_M3 input=WPG_Zug_M3', num2str(qq), ' interactive'];
end
system(cmd_str);
cd ..

%         2.1.2. catch the simulated forces along the loading time
step = 'Step-1';
req = 'FIXER_1-1,Node FIXER_1-1.1,RF2';
currentPath = cd;
Path = [currentPath,'\sim'];                                   % <===   get the forces   %

switch n 
    case 1 
        OdbFile = 'WPG_sim_M1.odb';                     %%%%%%\\----------------//%
    case 2 
        OdbFile='WPG_sim_M2.odb';                     %%%%%%\\----------------//%     
    case 3
        OdbFile='WPG_sim_M3.odb';                     %%%%%%\\----------------//%     
end

get_history_output(Path,OdbFile,step,req); 
ar = load('RF2.txt','-ascii');
f_fixer = ar(:,2)*1000;    % force in [mN]
tsim = ar(:,1);    % get the series of time
u_fixer = tsim *24;
e_fixer = u_fixer/40*100;    % strain in [%]

delete *.txt
cd ..
cd ..


end