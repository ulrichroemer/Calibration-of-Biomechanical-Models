function [posteriorMean,posteriorSample] = uqlab_update_Bayes(BayesOpts)

    myBayesianAnalysis = uq_createAnalysis(BayesOpts);
    uq_postProcessInversion(myBayesianAnalysis)
    posteriorSampleWalkers = myBayesianAnalysis.Results.PostProc.PostSample;
    chainSize = size(posteriorSampleWalkers,1)*size(posteriorSampleWalkers,3);
    numberParameters = length(BayesOpts.Prior.nonConst);
    posteriorSample = [];
    for i = 1:numberParameters
        parameterChain = reshape(squeeze(posteriorSampleWalkers(:,i,:)),[chainSize,1]);
        posteriorSample = [posteriorSample,parameterChain];
    end
    posteriorMean = myBayesianAnalysis.Results.PostProc.PointEstimate.X{1};

end