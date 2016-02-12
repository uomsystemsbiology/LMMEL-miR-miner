function [ arrayOutputGOTerms ] = matchHumanGOTermsToNums( arrayInputGONums, structFuncSettings )
%MATCHHUMANGOTERMSTONUMS Summary of this function goes here
%   Detailed explanation goes here

    [ arrayAllGONums, arrayAllGOTerms ] = extractHumanGOTermMappings( structFuncSettings );
    
    arrayOutputGOTerms = cell(length(arrayInputGONums),1);
    
    for iOutGOTerm = 1:length(arrayInputGONums),
        numInputGONum = arrayInputGONums(iOutGOTerm);
        numIndexforGOInFullArray = find(arrayAllGONums == numInputGONum);
        if isempty(numIndexforGOInFullArray),
            arrayOutputGOTerms{iOutGOTerm} = '-';
        else
            arrayOutputGOTerms{iOutGOTerm} = arrayAllGOTerms{numIndexforGOInFullArray};
        end
    end

end

