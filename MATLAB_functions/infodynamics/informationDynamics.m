function result = informationDynamics(vec1, vec2, whichFunction)
%% result = informationDynamics(vec1, vec2, whichFunction)
% This function is designed to take in two vectors of data and a string
%  specifying the function for calculating measures of pairwise statistical
%  assocation
%
%  Inputs:
%   - vec1: a vector of double precision values
%   - vec2: a vector of double precision values
%   - whichFunction: a string specifying the information dynamics metric,
%       expected values:
%       'TransferEntropyCalculatorKernel' - transfer entropy using the
%           kernel method within the JIDT
%       'TransferEntropyCalculatorKraskov' - transfer entropy using the
%           Kraskov estimator method within the JIDT
%       'MutualInfoCalculatorMultiVariateKernel' - mutual information using 
%           the kernel method within the JIDT
%       'MutualInfoCalculatorMultiVariateKraskov1' - mutual information
%           using the 1D Kraskov estimator method within the JIDT
%       'MutualInfoCalculatorMultiVariateKraskov2' - mutual information 
%           using the 2D Kraskov estimator method within the JIDT
%       'PearsCorr' - Pearson's correlation using the corr function within
%           the Statistics and Machine learning toolbox
%
%  Output:
%   - result: a double precision scalar specifying the calculated metric
%               
%  MATLAB Toolbox Dependencies:
%   - Statistics and Machine Learning Toolbox (for calc in 
%       informationDynamics)
%
%  Function dependencies:
%   - JIDT
% 
% This function was created by David Budden:
%   davidmarkbudden@gmail.com
%
% Last Updated: 03/03/16
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform pre-processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
 % -> Statistics and Machine Learning Toolbox (for corr)
%
% 
%TRANSFERENTROPY Summary of this function goes here
%   Detailed explanation goes here

%check that the input data are column vectors
if size(vec1,1) == 1 && size(vec1,2) > 1,
    disp('Warning; input vector one being used to calculate mutual information may need to be transposed');
end
if size(vec2,1) == 1 && size(vec2,2) > 1,
    disp('Warning; input vector two being used to calculate mutual information may need to be transposed');
end

if (size(vec1,2) == size(vec2,2)),
    numExtraDim = size(vec1,2);
end

%initialise the Java objects for the infoDynamics package
teCalc=initInfoDynamics(whichFunction,vec1,vec2);

if strcmp(whichFunction, 'TransferEntropyCalculatorKernel') || strcmp(whichFunction, 'TransferEntropyCalculatorKraskov')
    teCalc.initialise(); % Initialise leaving the parameters the same
    teCalc.setObservations(vec1, vec2);
    result = teCalc.computeAverageLocalOfObservations();
elseif strcmp(whichFunction, 'MutualInfoCalculatorMultiVariateKernel')
    teCalc.initialise(numExtraDim,numExtraDim);
    teCalc.setObservations(vec1, vec2);
    result = teCalc.computeAverageLocalOfObservations();
elseif strcmp(whichFunction, 'MutualInfoCalculatorMultiVariateKraskov1') || strcmp(whichFunction, 'MutualInfoCalculatorMultiVariateKraskov2')
    teCalc.initialise(numExtraDim,numExtraDim);
    teCalc.setObservations(vec1, vec2);
    result = teCalc.computeAverageLocalOfObservations();
elseif strcmp(whichFunction, 'PearsCorr')
    result = corr(vec1, vec2);
end

