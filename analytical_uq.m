% Iterative Bayesian update with model error and analytical function (Section 3.1)
% Reference: U. Roemer et al., "Surrogate-Based Bayesian Calibration of Biomechanical Models with Isotropic Material Behavior".

clearvars
uqlab

addpath(genpath(pwd))
close all
rng('default')
verbose = 0;
export = 0;

numberParameters = 2;
thetaTrue = [0.5,0.5];
noiseSD = 0.1;
xPointsPlot = linspace(0,1,200);

numberDataPoints = 10;    % Figure 5, left
%numberDataPoints = 30;    % Figure 5, right

xPoints = linspace(0,1,numberDataPoints);

noise = repmat(noiseSD,[1,numberDataPoints]);
data = analytical_model(thetaTrue,xPoints)' + noise.*randn(1,numberDataPoints);

if verbose
    plot(xPointsPlot,analytical_model([0.5,0.5],xPointsPlot),xPoints,data,'r*')
end

%% surrogate modeling

pceOpts.Degree = 10;       
numberSurrogateSamples = 1e5;

model = @(theta) analytical_model(theta,xPoints);

for i = 1:numberParameters
    IOpts.Marginals(i).Type = 'Uniform';
    IOpts.Marginals(i).Parameters = [0,1];
end

X_design = lhsdesign(numberSurrogateSamples,numberParameters);
Y_design = model(X_design);

number_KLE_modes = 5;
[eigenvalues,eigenvectorMatrix,eigenvalueMatrix] = uqlab_KLE(Y_design',number_KLE_modes);
xi_design = uqlab_KLE_sample(eigenvalues,eigenvectorMatrix,Y_design');    % project model data on KLE-space

pceOpts.ExpDesign.X = X_design;
pceOpts.ExpDesign.Y = xi_design;
pceOpts.Type = 'Metamodel';    
pceOpts.MetaType = 'PCE';            
pceOpts.Method = 'OLS';
pceOpts.Display = 'quiet';

myInput = uq_createInput(IOpts);
pce_metamodel = uq_createModel(pceOpts);  
surrogateModel = @(theta) uqlab_KLE_surrogate(eigenvectorMatrix,eigenvalueMatrix,mean(Y_design,2),uq_evalModel(pce_metamodel,theta));

%% standard Bayes (sigma fixed)

numberWalkers = 2*numberParameters;  
startValuesAIES = [unifrnd(0,1,numberWalkers,1), unifrnd(0,1,numberWalkers,1)];
samplesPerChain = 300;
burnInLength = 50;
stepSizeChain = 2; 

priorPDF = @(x) unifpdf(x(:,1),0,1).*unifpdf(x(:,2),0,1);
logPosterior = @(x) userLogLikelihood(x,noise,data,surrogateModel) + log(priorPDF(x));

samplesAIES = customAIES(startValuesAIES,logPosterior,samplesPerChain,stepSizeChain,burnInLength);

% histograms with true values and data vs posterior predictive 
if verbose
  
    figure 
    subplot(1,2,1)
    histogram(samplesAIES(:,1),'Normalization','Probability')
    hold on 
    plot(0.5,0,'r*')
    xlim([0,1])
    xlabel('parameter 1')

    subplot(1,2,2)
    histogram(samplesAIES(:,2),'Normalization','Probability')
    hold on 
    plot(0.5,0,'r*')
    xlim([0,1])
    xlabel('parameter 2')

    figure 
    posteriorPushForward = surrogateModel(samplesAIES);
    posteriorPredictive = posteriorPushForward' + repmat(noise,[length(samplesAIES),1]).*randn(length(samplesAIES),numberDataPoints);
    [posteriorPredictiveLower,posteriorPredictiveUpper] = getQuantiles(posteriorPredictive,0.99);
    plot(xPoints,data,'r*',xPoints,surrogateModel(mean(samplesAIES)),'b-',...
          xPoints,posteriorPredictiveLower,'k--',xPoints,posteriorPredictiveUpper,'k--')
    hold on 
    fill([xPoints, fliplr(xPoints)],[posteriorPredictiveLower, fliplr(posteriorPredictiveUpper)],'b', 'facealpha',0.05)
  
end

%% standard Bayes (sigma estimated)

numberWalkers = 2*(numberParameters+1);  
startValuesAIES = [unifrnd(0,1,numberWalkers,1), unifrnd(0,1,numberWalkers,1), unifrnd(0,2,numberWalkers,1)];
samplesPerChain = 200;
burnInLength = 40;
stepSizeChain = 2; 

logLikelihood = @(theta) userLogLikelihoodDiscrepancy(theta,data,surrogateModel);

priorPDF = @(x) unifpdf(x(:,1),0,1).*unifpdf(x(:,2),0,1).*unifpdf(x(:,3),0,2);
logPosterior = @(x) logLikelihood(x) + log(priorPDF(x));

samplesAIESSigma = customAIES(startValuesAIES,logPosterior,samplesPerChain,stepSizeChain,burnInLength);

%% iterative Bayesian with surrogate and residual model error 

% redefine prior and likelihood: include an additional parameter for remaining model error 
numberWalkers = 2*(numberParameters+1);  
startValuesAIES = [unifrnd(0,1,numberWalkers,1), unifrnd(0,1,numberWalkers,1),unifrnd(0,1,numberWalkers,1)];
samplesPerChain = 200;
burnInLength = 40;
stepSizeChain = 2; 
priorPDF = @(x) 0.5*unifpdf(x(:,1),0,1).*unifpdf(x(:,2),0,1).*unifpdf(x(:,3),0,1);

numberIterations = 10;
posteriorSampleIterations = cell(numberIterations,1);
modelError = zeros(numberDataPoints,length(samplesAIES));

for i = 1:numberIterations
    
    augmentedLogLikelihood = @(theta) userLogLikelihoodIterative(theta,modelError,data,surrogateModel);
    
    logPosterior = @(theta) augmentedLogLikelihood(theta) + log(priorPDF(theta));
    augmentedPosteriorSample = customAIES(startValuesAIES,logPosterior,samplesPerChain,stepSizeChain,burnInLength);
    posteriorSample = augmentedPosteriorSample(:,1:numberParameters);
    
    posteriorSampleIterations{i} = posteriorSample;
    modelResponseAtPosterior = model(posteriorSample);
    modelError = (modelResponseAtPosterior - surrogateModel(posteriorSample));
    if i>1
        disp(['mean theta_1: (i) ',num2str(mean(posteriorSampleIterations{i}(:,1))),' (i-1) ',num2str(mean(posteriorSampleIterations{i-1}(:,1)))])
        disp(['mean theta_2: (i) ',num2str(mean(posteriorSampleIterations{i}(:,2))),' (i-1) ',num2str(mean(posteriorSampleIterations{i-1}(:,2)))])
        %disp(['model error: mean ',num2str(norm(mean(modelError,2))),' std ',num2str(norm(std(modelError,2)))])
    end
    
end

samplesAIESIterative = posteriorSample;
mean(samplesAIESIterative)

%% plotting

% plot estimated parameters with credible intervals
if verbose
    figure
    plot(thetaTrue(1),thetaTrue(2),'r*','MarkerSize',5)
    
    posteriorMean = mean(samplesAIES,1);
    posteriorStd = std(samplesAIES,1);
    ellipse(posteriorMean(1),posteriorMean(2),posteriorStd(1),posteriorStd(2),'b-')
    

    posteriorMean = mean(samplesAIESSigma,1);
    posteriorStd = std(samplesAIESSigma,1);
    ellipse(posteriorMean(1),posteriorMean(2),posteriorStd(1),posteriorStd(2),'b--')
    
    posteriorMean = mean(samplesAIESIterative,1);
    posteriorStd = std(samplesAIESIterative,1);
    ellipse(posteriorMean(1),posteriorMean(2),posteriorStd(1),posteriorStd(2),'b-.')
    
    legend('true value','standard Bayes','standard Bayes (est. sigma)','iterative Bayes (est. sigma)','Location','southeast')
    xlim([-0.2,1])
    ylim([-0.2,1])
    
    set(gcf,'units','points','position',[200,200,470*0.5,470*0.5/1.618])
    set(gca,'FontSize',10)
end
    

if export 
    savefig(['analytical_',num2str(numberDataPoints),'.fig'])
    print(['analytical_parameters_iteration_',num2str(numberDataPoints)],'-dpng','-r400')
end

% plot surrogate error against data error
if verbose
    sigmaE = mean(samplesAIESSigma(:,3));

    figure
    plot(xPoints,sigmaE*ones(size(xPoints)),'r-',xPoints,noise,'r--',xPoints,abs(mean(modelError,2)),'b-')
    
    figure
    plot(xPoints,data,'r-',xPoints,surrogateModel(mean(samplesAIESSigma(:,1:2))))
end

% plot posterior vs data
if verbose 

    figure 

    % standard
    samplesPlot = samplesAIES;
    posteriorPredictive = surrogateModel(samplesPlot)';
    [posteriorPredictiveLower,posteriorPredictiveUpper] = getQuantiles(posteriorPredictive,0.95);
    plot(xPoints,data,'r*',xPoints,surrogateModel(mean(samplesPlot,1))','b-',...
          xPoints,posteriorPredictiveLower,'k--',xPoints,posteriorPredictiveUpper,'k--')
    hold on 
    fill([xPoints, fliplr(xPoints)],[posteriorPredictiveLower, fliplr(posteriorPredictiveUpper)],'b', 'facealpha',0.05)

    % corrected
    samplesPlot = samplesAIESIterative; 
    posteriorPredictive = surrogateModel(samplesPlot)';
    [posteriorPredictiveLower,posteriorPredictiveUpper] = getQuantiles(posteriorPredictive,0.95);
    plot(xPoints,data - mean(modelError,2)','g*',xPoints,surrogateModel(mean(samplesPlot,1))','g-',...
          xPoints,posteriorPredictiveLower,'g--',xPoints,posteriorPredictiveUpper,'g--')
    hold on 
    fill([xPoints, fliplr(xPoints)],[posteriorPredictiveLower, fliplr(posteriorPredictiveUpper)],'g', 'facealpha',0.05)
end

%% aux functions
function h = ellipse(x,y,rx,ry,string)
    hold on
    th = 0:pi/50:2*pi;
    xunit = rx * cos(th) + x;
    yunit = ry * sin(th) + y;
    h = plot(xunit, yunit, string);
    hold off
end

function L = userLogLikelihood(theta,noise,data,surrogateModel)
   
    covariance = diag(noise).^2;
    L = -(1/2)*log(det(covariance)) -(1/2)*sum((surrogateModel(theta(1:2))'-data).*(covariance\(surrogateModel(theta(1:2))'-data)')',2);  

end

function L = userLogLikelihoodDiscrepancy(theta,data,surrogateModel)
   
    covariance = diag((theta(3)).^2*ones(size(data)));
    L = -(1/2)*log(det(covariance)) -(1/2)*sum((surrogateModel(theta(1:2))'-data).*(covariance\(surrogateModel(theta(1:2))'-data)')',2);  

end

function L = userLogLikelihoodIterative(theta,modelError,data,surrogateModel)
    
    augmentedData = (data - mean(modelError,2)');
    augmentedCovariance = diag((theta(3)).^2*ones(size(data))) + cov(modelError');
    L = -(1/2)*log(det(augmentedCovariance)) -(1/2)*sum((surrogateModel(theta(1:2))'-augmentedData).*(augmentedCovariance\(surrogateModel(theta(1:2))'-augmentedData)')',2);  

end

