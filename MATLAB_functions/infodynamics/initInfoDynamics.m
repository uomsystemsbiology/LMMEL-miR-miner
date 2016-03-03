function teCalc = initInfoDynamics(whichFunction, vec1, vec2)
%DAVIDINITIALISETE Summary of this function goes here
%   Detailed explanation goes here

if (size(vec1,2) == size(vec2,2)),
    numExtraDim = size(vec1,2);
end

if strcmp(whichFunction, 'TransferEntropyCalculatorKernel')
    teCalc=javaObject('infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel');
    teCalc.setProperty('NORMALISE_PROP_NAME', 'false'); % Normalise the individual variables
    teCalc.initialise(1, 0.5); % Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
elseif strcmp(whichFunction, 'TransferEntropyCalculatorKraskov')
    teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');    
    teCalc.setProperty('k', '4'); % Use Kraskov parameter K=4 for 4 nearest points
    teCalc.initialise(1); % Use history length 1 (Schreiber k=1)
elseif strcmp(whichFunction, 'MutualInfoCalculatorMultiVariateKernel')
    teCalc=javaObject('infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel');
    teCalc.setProperty('NORMALISE_PROP_NAME', 'false'); % Normalise the individual variables
    teCalc.initialise(numExtraDim,numExtraDim);
elseif strcmp(whichFunction, 'MutualInfoCalculatorMultiVariateKraskov1');
    teCalc=javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1');
    teCalc.initialise(numExtraDim,numExtraDim);
elseif strcmp(whichFunction, 'MutualInfoCalculatorMultiVariateKraskov2');
    teCalc=javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');
    teCalc.initialise(numExtraDim,numExtraDim);
else
    teCalc = -1;
end

end

