function [samplesEMCMC] = customAIES(startValues,logPosterior,samplesPerChain,stepSizeChain,burnInLength)

  outputEMCMC = EMCMCsampler(startValues,logPosterior,samplesPerChain,'StepSize',stepSizeChain,'burnin',burnInLength);

  % re-format the output sample
  samplesEMCMC = outputEMCMC.samples;
  samplesEMCMC = permute(samplesEMCMC, [2 1 3]);
  samplesEMCMC = samplesEMCMC(:,:)';

end