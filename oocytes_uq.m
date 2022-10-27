% Iterative Bayesian update with model error and gel example (Section 3.3)
% Reference: U. Roemer et al., "Surrogate-Based Bayesian Calibration of Biomechanical Models with Isotropic Material Behavior".

clearvars

addpath(genpath(pwd))
close all
rng('default')
verbose = 1;
export = 0;
uqlab

numberParameters = 4;    % number of parameters (code only works for 4 parameters)
xPoints = 1:49;    % loading indices
id = '5dc3066b-1604-40d6-97ef-401b1b3ae3b6';    % 500 design points (4 parameters)
F_ExperimentalDesign = csvread(['data/exp_design_output_material_id_',id,'.csv']);
F_ExperimentalDesign = F_ExperimentalDesign(:,xPoints)';
G_ExperimentalDesign = log(F_ExperimentalDesign)';

theta_ExperimentalDesign = csvread(['data/exp_design_input_material_id_',id,'.csv']);
theta_ExperimentalDesign = rescale(theta_ExperimentalDesign);

[data,noise] = get_data();
numberDataPoints = length(data);

%% setup surrogate

number_KLE_modes = 3;
polynomialDegree = 3;

[eigenvalues,eigenvectorMatrix,eigenvalueMatrix] = uqlab_KLE(G_ExperimentalDesign,number_KLE_modes);
xi_design = uqlab_KLE_sample(eigenvalues,eigenvectorMatrix,G_ExperimentalDesign);   

pceOpts.ExpDesign.X = theta_ExperimentalDesign;
pceOpts.ExpDesign.Y = xi_design;
pceOpts.Type = 'Metamodel';    
pceOpts.MetaType = 'PCE';            
pceOpts.Method = 'OLS';

for i = 1:numberParameters
    IOpts.Marginals(i).Type = 'Uniform';
    IOpts.Marginals(i).Parameters = [-1,1];    
end
myInput = uq_createInput(IOpts);
pceOpts.Degree = polynomialDegree;       
pce_metamodel = uq_createModel(pceOpts);   
surrogateModel = @(theta) exp(uqlab_KLE_surrogate(eigenvectorMatrix,eigenvalueMatrix,mean(G_ExperimentalDesign,1)',uq_evalModel(pce_metamodel,theta)));

% reference model 
number_KLE_modesRef = 4;
polynomialDegreeRef = 4;
pceOptsRef = pceOpts;
pceOptsRef.Degree = polynomialDegreeRef;  

[eigenvaluesRef,eigenvectorMatrixRef,eigenvalueMatrixRef] = uqlab_KLE(G_ExperimentalDesign,number_KLE_modesRef);
xi_designRef = uqlab_KLE_sample(eigenvaluesRef,eigenvectorMatrixRef,G_ExperimentalDesign);   
pceOptsRef.ExpDesign.Y = xi_designRef;
pce_refMetamodel = uq_createModel(pceOptsRef);   
referenceModel = @(theta) exp(uqlab_KLE_surrogate(eigenvectorMatrixRef,eigenvalueMatrixRef,mean(G_ExperimentalDesign,1)',uq_evalModel(pce_refMetamodel,theta)));

%% standard Bayesian update

numberWalkers = 4*numberParameters;  
startValuesAIES = [unifrnd(-1,1,numberWalkers,1), unifrnd(-1,1,numberWalkers,1), unifrnd(-1,1,numberWalkers,1), unifrnd(-1,1,numberWalkers,1)];
samplesPerChain = 250;
burnInLength = 30;
stepSizeChain = 1; 

priorPDF = @(x) unifpdf(x(:,1),-1,1).*unifpdf(x(:,2),-1,1).*unifpdf(x(:,3),-1,1).*unifpdf(x(:,4),-1,1);
logPosterior = @(x) userLogLikelihood(x,noise,data,surrogateModel) + log(priorPDF(x));

samplesAIES = customAIES(startValuesAIES,logPosterior,samplesPerChain,stepSizeChain,burnInLength);

%% iterative Bayesian with surrogate and residual model error

% redefine prior and likelihood: include an additional parameter for remaining model error 
numberWalkers = 2*(numberParameters + 1);  
startValuesAIES = [unifrnd(-1,1,numberWalkers,1), unifrnd(-1,1,numberWalkers,1),unifrnd(-1,1,numberWalkers,1),unifrnd(-1,1,numberWalkers,1),unifrnd(0,1,numberWalkers,1)];
samplesPerChain = 150;
burnInLength = 30;
stepSizeChain = 1; 
priorPDF = @(x) unifpdf(x(:,1),-1,1).*unifpdf(x(:,2),-1,1).*unifpdf(x(:,3),-1,1).*unifpdf(x(:,4),-1,1).*unifpdf(x(:,5),0,1);

numberIterations = 10;
posteriorSampleIterations = cell(numberIterations,1);
modelError = zeros(numberDataPoints,length(samplesAIES));

for i = 1:numberIterations
    
    augmentedLogLikelihood = @(theta) userLogLikelihoodIterative(theta,modelError,noise,data,surrogateModel);
    
    logPosteriorAugmented = @(theta) augmentedLogLikelihood(theta) + log(priorPDF(theta));
    augmentedPosteriorSample = customAIES(startValuesAIES,logPosteriorAugmented,samplesPerChain,stepSizeChain,burnInLength);
    posteriorSample = augmentedPosteriorSample(:,1:numberParameters);
    
    if(min(data - mean(modelError,2)')<0)
        disp('surrogate too coarse')
    end
    posteriorSampleIterations{i} = posteriorSample;
    modelResponseAtPosterior = referenceModel(posteriorSample);
    modelError = (modelResponseAtPosterior - surrogateModel(posteriorSample));
    if i>1
        disp(['mean theta_1: (i) ',num2str(mean(posteriorSampleIterations{i}(:,1))),' (i-1) ',num2str(mean(posteriorSampleIterations{i-1}(:,1)))])
        disp(['mean theta_2: (i) ',num2str(mean(posteriorSampleIterations{i}(:,2))),' (i-1) ',num2str(mean(posteriorSampleIterations{i-1}(:,2)))])
        disp(['mean theta_3: (i) ',num2str(mean(posteriorSampleIterations{i}(:,3))),' (i-1) ',num2str(mean(posteriorSampleIterations{i-1}(:,3)))])
        disp(['mean theta_4: (i) ',num2str(mean(posteriorSampleIterations{i}(:,4))),' (i-1) ',num2str(mean(posteriorSampleIterations{i-1}(:,4)))])
        disp('')
    end
    
end

samplesAIESIterative = posteriorSample;

%% plot posterior vs data

if verbose 

    figure 

    % standard
    samplesPlot = samplesAIES;
    [~,indMAP] = max(logPosterior(samplesPlot));
    posteriorPushForward = surrogateModel(samplesPlot);
    posteriorPredictive = posteriorPushForward';
    [posteriorPredictiveLower,posteriorPredictiveUpper] = getQuantiles(posteriorPredictive,0.95);
    plot(xPoints,data,'r*',xPoints,surrogateModel(samplesPlot(indMAP,:)),'b-',...
          xPoints,posteriorPredictiveLower,'k--',xPoints,posteriorPredictiveUpper,'k--')
    hold on 
    fill([xPoints, fliplr(xPoints)],[posteriorPredictiveLower, fliplr(posteriorPredictiveUpper)],'b', 'facealpha',0.05)


    figure 
    samplesPlot = samplesAIESIterative;  
    [~,indMAP] = max(logPosteriorAugmented([samplesPlot,zeros(size(samplesPlot,1),1)]));
    posteriorPushForward = surrogateModel(samplesPlot);
    posteriorPredictive = posteriorPushForward';
    [posteriorPredictiveLower,posteriorPredictiveUpper] = getQuantiles(posteriorPredictive,0.95);
    plot(xPoints,data - mean(modelError,2)','g*',xPoints,surrogateModel(samplesPlot(indMAP,:)),'g-',...
          xPoints,posteriorPredictiveLower,'g--',xPoints,posteriorPredictiveUpper,'g--')
    hold on 
    fill([xPoints, fliplr(xPoints)],[posteriorPredictiveLower, fliplr(posteriorPredictiveUpper)],'g', 'facealpha',0.05)
end

%% parameter statistics
CI_posteriorSample = 10.^(posteriorSample);
C_quantiles = plims(CI_posteriorSample,[0.05,0.5,0.95]);    

%% aux functions

function L = userLogLikelihood(theta,noise,data,surrogateModel)
   
    covariance = diag(noise).^2;
    L = -(1/2)*log(det(covariance)) -(1/2)*sum((surrogateModel(theta)'-data).*(covariance\(surrogateModel(theta)'-data)')',2);  

end

function L = userLogLikelihoodIterative(theta,modelError,noise,data,surrogateModel)
    
    augmentedData = (data - mean(modelError,2)');
    augmentedCovariance = diag((theta(5)).^2*ones(size(data))) + cov(modelError');
    L = -(1/2)*log(det(augmentedCovariance)) - (1/2)*sum((surrogateModel(theta(1:4))'-augmentedData).*(augmentedCovariance\(surrogateModel(theta(1:4))'-augmentedData)')',2);  

end


function [theta_rs] = rescale(theta)

    % lin-scale
    mu_min = 1e-3;
    mu_max = 30e-3;

    alpha_min = 5;
    alpha_max = 25;
    
    M = size(theta,2);
    
    if M == 2
        theta_rs_1 = 2*((theta(:,1) - mu_min)/(mu_max - mu_min) - 1) + 1;
        theta_rs_2 = 2*((theta(:,2) - alpha_min)/(alpha_max - alpha_min) - 1) + 1;
        theta_rs = [theta_rs_1, theta_rs_2];
    end
    
    if M == 4
        theta_rs_1 = 2*((theta(:,1) - mu_min)/(mu_max - mu_min) - 1) + 1;
        theta_rs_2 = 2*((theta(:,2) - alpha_min)/(alpha_max - alpha_min) - 1) + 1;
        theta_rs_3 = 2*((theta(:,3) - mu_min)/(mu_max - mu_min) - 1) + 1;
        theta_rs_4 = 2*((theta(:,4) - alpha_min)/(alpha_max - alpha_min) - 1) + 1;
        theta_rs = [theta_rs_1, theta_rs_2, theta_rs_3, theta_rs_4];
    end
    
   
end


function [r] = upscale(s)

    M = size(s,2);

    % lin-scale
    mu_min = 1e-3;
    mu_max = 30e-3;

    alpha_min = 5;
    alpha_max = 25;
    
    if M == 2
        r(:,1) = mu_min + (mu_max-mu_min)*(1+s(:,1))/2;
        r(:,2) = alpha_min + (alpha_max-alpha_min)*(1+s(:,2))/2;
    end
    
    if M == 4
        r(:,1) = mu_min + (mu_max-mu_min)*(1+s(:,1))/2;
        r(:,2) = alpha_min + (alpha_max-alpha_min)*(1+s(:,2))/2;
        
        r(:,3) = mu_min + (mu_max-mu_min)*(1+s(:,3))/2;
        r(:,4) = alpha_min + (alpha_max-alpha_min)*(1+s(:,4))/2;    
    end
   
end


function [b,sig] = get_data()
    
    data = load('data/daten_fig5_indent_raw.txt');
    b = data(:,3);
    blow = data(:,2);
    bup = data(:,4);
    eps = data(:,1);
    
    x_int = linspace(2,58,50);
    
    b = abs(b');
    sig = abs(bup-blow)/2;
    sig = sig';
    
    % change from eps to loading increment 
    b = interp1(-100*eps,b,x_int);
    sig = interp1(-100*eps,sig,x_int);
    
    b = b(2:end);
    sig = sig(2:end);

end