% Iterative Bayesian update with model error and gel example (Section 3.2)
% Reference: U. Roemer et al., "Surrogate-Based Bayesian Calibration of Biomechanical Models with Isotropic Material Behavior".

clearvars

addpath(genpath(pwd))
close all
rng('default')
verbose = 1;
export = 0;
uqlab

numberParameters = 3;    % number of parameters (code is set up for 3 parameters)
xPoints = [3.0000, 6.0000, 9.0000, 12.0000, 15.0000, 18.0000, 21.0000, 24.0000, 27.0000, 30.0000, 33.0000, 36.0000, 39.0000, 42.0000, 45.0000, 48.0000, 51.0000, 54.0000, 57.0000, 60.0000];

[xMeas,fMeas,sigma_meas] = get_data();
data = interp1(xMeas,fMeas,xPoints);
noise = sigma_meas(2:end)';
numberDataPoints = length(data);
xPointsData = xMeas(2:end)';

%% setup surrogate

id = '0c3f341d-2682-493c-9674-43cc5c0c6b4e';    % id for data with 1000 sample points
id_in = ['data/gel_exp_design_input_material_id_',id,'.csv'];
id_out = ['data/gel_simulated_data_id_',id,'.csv'];

F_ExperimentalDesign = csvread(id_out);
F_ExperimentalDesign = F_ExperimentalDesign(2:end,2:end)';    % first column contains displacements -> remove
G_ExperimentalDesign = log(F_ExperimentalDesign);

theta_ExperimentalDesign = csvread(id_in);
theta_ExperimentalDesign = rescale(theta_ExperimentalDesign);

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
    IOpts.Marginals(i).Parameters = [-5,1];    
end
myInput = uq_createInput(IOpts);
pceOpts.Degree = polynomialDegree;       
pce_metamodel = uq_createModel(pceOpts);   
surrogateModel = @(theta) exp(uqlab_KLE_surrogate(eigenvectorMatrix,eigenvalueMatrix,mean(G_ExperimentalDesign,1)',uq_evalModel(pce_metamodel,theta)));

% reference model 
number_KLE_modesRef = 8;
polynomialDegreeRef = 9;
pceOptsRef = pceOpts;
pceOptsRef.Degree = polynomialDegreeRef;  

[eigenvaluesRef,eigenvectorMatrixRef,eigenvalueMatrixRef] = uqlab_KLE(G_ExperimentalDesign,number_KLE_modesRef);
xi_designRef = uqlab_KLE_sample(eigenvaluesRef,eigenvectorMatrixRef,G_ExperimentalDesign);   
pceOptsRef.ExpDesign.Y = xi_designRef;
pce_refMetamodel = uq_createModel(pceOptsRef);   
referenceModel = @(theta) exp(uqlab_KLE_surrogate(eigenvectorMatrixRef,eigenvalueMatrixRef,mean(G_ExperimentalDesign,1)',uq_evalModel(pce_refMetamodel,theta)));

%% standard Bayesian update

numberWalkers = 2*numberParameters;  
startValuesAIES = [unifrnd(-5,1,numberWalkers,1), unifrnd(-5,1,numberWalkers,1), unifrnd(-5,1,numberWalkers,1)];
samplesPerChain = 150;
burnInLength = 30;
stepSizeChain = 2; 

priorPDF = @(x) unifpdf(x(:,1),-5,1).*unifpdf(x(:,2),-5,1).*unifpdf(x(:,3),-5,1);
logPosterior = @(x) userLogLikelihood(x,noise,data,surrogateModel) + log(priorPDF(x));

samplesAIES = customAIES(startValuesAIES,logPosterior,samplesPerChain,stepSizeChain,burnInLength);
mean(samplesAIES)

%% iterative Bayesian with surrogate and residual model error

% redefine prior and likelihood: include an additional parameter for remaining model error 
numberWalkers = 2*(numberParameters);  
startValuesAIES = [unifrnd(-5,1,numberWalkers,1), unifrnd(-5,1,numberWalkers,1),unifrnd(-5,1,numberWalkers,1)];
samplesPerChain = 150;
burnInLength = 30;
stepSizeChain = 2; 
priorPDF = @(x) unifpdf(x(:,1),-5,1).*unifpdf(x(:,2),-5,1).*unifpdf(x(:,3),-5,1);

numberIterations = 5;
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
        disp('')
        %disp(['model error: mean ',num2str(norm(mean(modelError,2))),' std ',num2str(norm(std(modelError,2)))])
    end
    
end

samplesAIESIterative = posteriorSample;
mean(samplesAIESIterative)

%% plot posterior vs data
close all

if verbose 
    figure 

    % standard
    samplesPlot = samplesAIES;
    posteriorPushForward = surrogateModel(samplesPlot);
    posteriorPredictive = posteriorPushForward';
    [posteriorPredictiveLower,posteriorPredictiveUpper] = getQuantiles(posteriorPredictive,0.95);
    BayesEstimate = surrogateModel(mean(samplesPlot,1));
    plot(xPointsData,data,'r*')    

    % add iterative component
    samplesPlot = samplesAIESIterative;  
    posteriorPushForward = surrogateModel(samplesPlot);
    posteriorPredictive = posteriorPushForward';
    [posteriorPredictiveLower,posteriorPredictiveUpper] = getQuantiles(posteriorPredictive,0.95)
    BayesEstimate = surrogateModel(mean(samplesPlot,1));
    plot(xPointsData,data,'r*',xPointsData,data - mean(modelError,2)','g*',xPointsData,BayesEstimate + mean(modelError,2),'g-',...
          xPointsData,posteriorPredictiveLower + mean(modelError,2)','g--',xPointsData,posteriorPredictiveUpper + mean(modelError,2)','g--')
    hold on 
    fill([xPointsData, fliplr(xPointsData)],[posteriorPredictiveLower + mean(modelError,2)', fliplr(posteriorPredictiveUpper + mean(modelError,2)')],'g', 'facealpha',0.05)
end

if export 
    h = legend('data','posterior mean','90% credible intervals'); 
    set(gcf,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
    set(h, 'Position', [0.4, 0.7, 0.05, .1])
    print('Gel_posterior','-dpng','-r400')
end

%% plot model error 

if verbose 
    priorSample = theta_ExperimentalDesign;
    priorError = abs(mean(referenceModel(priorSample) - surrogateModel(priorSample),2));
    posteriorError = abs(mean(modelError,2));

    priorErrorStored33 = 1e3*[0.0060    0.0920    0.2007    0.3002    0.3954    0.4954    0.6084    0.7414    0.8997    1.0876    1.3085 1.5653    1.8609    2.1982    2.5812    3.0087    3.4862    4.0150    4.5968    5.2323];
    posteriorErrorStored33 = [0.1293    0.3776    1.1280    1.9221    2.7370    3.5383    4.2831    4.9366    5.4820    5.9163    6.2462 6.4851    6.6479    6.7497    6.8053    6.8383    6.8569    6.8759    6.9077    6.9628];

    posteriorErrorStored55 = [0.1557    0.1564    0.0438    0.3726    0.7599    1.1595    1.5458    1.9074    2.2422    2.5525    2.8430 3.1183    3.3837    3.6452    3.9092    4.1703    4.4391    4.7169    5.0058    5.3076];
    priorErrorStored55 = [2.1361    3.9238    2.4851    7.4923   17.4736   29.0707   40.8574   52.9053   65.8590   80.4811   97.4768 117.4128  140.7565  167.8000  198.5717  232.0678  267.8186  304.3434  339.6926  371.3430];

    posteriorErrorStored77 = [0.0624    0.1311    0.1038    0.0205    0.0930    0.2161    0.3357    0.4451    0.5423    0.6286    0.7065 0.7789    0.8489    0.9195    0.9944    1.0713    1.1557    1.2484    1.3507    1.4636];
    priorErrorStored77 = [0.4506    1.1002    1.8052    1.7402    2.2194    3.3120    4.8403    6.5922    8.4701   10.4573   12.6025 15.0416   17.9899   21.6473   26.0799   31.4119   37.3584   43.4461   48.8452   52.2765];


    figure 
    semilogy(xPointsData,posteriorErrorStored33,'r-',xPointsData,posteriorErrorStored55,'r--',xPointsData,posteriorErrorStored77,'r.-')
end

% export
export = 1;
if export 
    set(gcf,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
    h = legend({'$T=3,P=3$','$T=5,P=5$','$T=7,P=7$'},'Interpreter','latex');
    set(h, 'Position', [0.7, 0.2, 0.05, .1])
    %print('Gel_prior_error','-dpng','-r400')
    print('Gel_posterior_error','-dpng','-r400')
end

%% parameter statistics
CI_posteriorSample = 10.^(posteriorSample);
C_quantiles = plims(CI_posteriorSample,[0.05,0.5,0.95]);    

%% Sobol sensitivity analysis
rng('default')

ind = 3;
samplesMC = 1e5;
X = -5 + 6*rand(samplesMC,numberParameters);

% generate second sample, where part of the parameter set is resampled
Xu = X;
Xu(:,setdiff(1:numberParameters,ind)) = rand(samplesMC,numberParameters-1);

n = 20;
pceOptsI = pceOpts;    
pceOptsI.degree = 5;
pceOptsI.ExpDesign.X = theta_ExperimentalDesign;
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.Order = 1;

% with estimator
Y = surrogateModel(X);
Yu = surrogateModel(Xu);

genSobolNum = 0;
genSobolDenum = 0;
for i = 1:numberDataPoints
    
    genSobolNum = genSobolNum + mean(Y(i,:).*Yu(i,:)) - mean(Y(i,:))*mean(Yu(i,:));
    genSobolDenum = genSobolDenum + var(Y(i,:));
  
end

S_ind = genSobolNum/genSobolDenum


%% aux functions

function L = userLogLikelihood(theta,noise,data,surrogateModel)
   
    covariance = diag(noise).^2;
    L = -(1/2)*log(det(covariance)) -(1/2)*sum((surrogateModel(theta)'-data).*(covariance\(surrogateModel(theta)'-data)')',2);  

end

function L = userLogLikelihoodIterative(theta,modelError,noise,data,surrogateModel)
    
    augmentedData = (data - mean(modelError,2)');
    augmentedCovariance = diag(noise.^2) + cov(modelError');
    L = -(1/2)*log(det(augmentedCovariance)) - (1/2)*sum((surrogateModel(theta(1:3))'-augmentedData).*(augmentedCovariance\(surrogateModel(theta(1:3))'-augmentedData)')',2);  

end


function [theta_rs] = rescale(theta)

    theta_rs = log10(theta);
    
end

function [theta] = upscale(theta_rs)

    theta = 10.^(theta_rs);    
end

function [xMeas,fMeas,sigma_meas] = get_data()

    strainForce = xlsread('data/exp_WPG.xlsx');
    xMeas = strainForce(:,1);
    fMeas = mean(strainForce(:,2:end),2);
    sigma_meas = std(strainForce(:,2:end)')';
    xMeas = xMeas(1:6:end);
    fMeas = fMeas(1:6:end);
    sigma_meas = sigma_meas(1:6:end);
    
end

