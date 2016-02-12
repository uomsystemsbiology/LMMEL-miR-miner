function [ arrayOutputEnsGeneIDs, arrayOutputGOTerms, arrayOutputGONums ] = biomartEnsemblGeneToGeneOntology( arrayInputGOTerms, structFuncSettings )
%EXTRACTENTREZGENEINFO Summary of this function goes here
%   Detailed explanation goes here


    numGOTermsOfInterest = length(arrayInputGOTerms);
    
    stringEnsembleIDMatcher = 'ENSG';

    if structFuncSettings.performNewExtraction,

        disp([ char(9) '... extracting gene information from the Ensembl database']);
        
        %open the input file 
        fileInputData = fopen([ structFuncSettings.biomartFolder 'ensemblGene_to_GO.csv' ],'rt');
        arrayInputCell = textscan(fileInputData,'%s','delimiter','\n');
        fclose(fileInputData);
        arrayBioMartMapping = arrayInputCell{1};
        clear arrayInputCell;

        numRows = length(arrayBioMartMapping);
        
        arrayOutputRowIndices = [2:numRows];
        numOutputRows = length(arrayOutputRowIndices);
        
        arrayDatabaseEnsemblGeneIDs = zeros(numOutputRows,1,'uint32');
        arrayDatabaseGOTerms = cell(numOutputRows,1);
        arrayDatabaseGONums = zeros(numOutputRows,1, 'uint32');
        arrayGOOfIntFlag = false(numOutputRows,1);
        for iDatabaseRow = 1:numOutputRows,
            numRow = arrayOutputRowIndices(iDatabaseRow);
            stringRow = arrayBioMartMapping{numRow};
            arrayCommaPos = strfind(stringRow, ',');
            
            %if the GO term contains a comma, then it will be encased
            %within quotation marks
            arrayQuotationPos = strfind(stringRow, '"');
            if ~isempty(arrayQuotationPos),
                if length(arrayQuotationPos) == 2,
                    arrayCommaInQuotes = (arrayCommaPos > arrayQuotationPos(1)) & (arrayCommaPos < arrayQuotationPos(2));
                    arrayCommaPos = arrayCommaPos(~arrayCommaInQuotes);
                else
                    disp('warning: there are more than two quotation marks, the code needs to be fixed accordingly');
                end
            end
            
            %extract the Ensembl gene ID from the firstcolumn
            stringEnsemblGeneID = stringRow(1:(arrayCommaPos(1)-1));
            if ~isempty(stringEnsemblGeneID),
                numEnsemblGeneIDDouble = str2double(stringEnsemblGeneID((length(stringEnsembleIDMatcher)+1):end));
                if numEnsemblGeneIDDouble < (2^32 - 1),
                    arrayDatabaseEnsemblGeneIDs(iDatabaseRow) = uint32(numEnsemblGeneIDDouble);
                else
                    disp('warning: ensembl gene IDs exceed values that can be stored as uint32, code should be rewritten with uint64');
                end
            end
                
            %extract the GO term from the second column
            stringGOTerm = stringRow((arrayCommaPos(1)+1):(arrayCommaPos(2)-1));
            if ~isempty(stringGOTerm),
                %remove the quotation marks if they exist
                if strncmp(stringGOTerm(1), '"', 1),
                    stringGOTerm = stringGOTerm(2:end);
                end
                if strncmp(stringGOTerm(end), '"', 1),
                    stringGOTerm = stringGOTerm(1:(end-1));
                end
                arrayDatabaseGOTerms{iDatabaseRow} = stringGOTerm;
                
                %search for the GO terms of interest and set the flag
                for iGOSearchTerm = 1:numGOTermsOfInterest,
                    
                    arrayTermPosition = strfind(stringGOTerm, arrayInputGOTerms{iGOSearchTerm});
                    
                    if ~isempty(arrayTermPosition),
                        arrayGOOfIntFlag(iDatabaseRow) = true;
                    end
                    
                end
                
            else
                arrayDatabaseGOTerms{iDatabaseRow} = '-';
            end
            
            
            %extract the GO term index number from the third column
            stringGONum = stringRow((arrayCommaPos(2)+1):end);
            if ~isempty(stringGONum),
                
                if strncmp(stringGONum, 'GO:', 3),
                
                    numGODouble = str2double(stringGONum(4:end));
                    if numGODouble < (2^32 - 1),
                        arrayDatabaseGONums(iDatabaseRow) = uint32(numGODouble);
                    else
                        disp('warning: GO term numbers exceed values that can be stored as uint32, code should be rewritten with uint64');
                    end
                
                else
                    disp('warning: GO term number not in the expected format');
                end
                
            else
                %vector initialised as zero, do nothing
            end
                        
        end
                
        save([ structFuncSettings.biomartFolder 'biomartEnsemblGeneToGeneOntology.mat'], 'arrayDatabaseEnsemblGeneIDs', 'arrayDatabaseGOTerms', 'arrayDatabaseGONums', 'arrayGOOfIntFlag');
        
    else
        
        structBioMartData = load([ structFuncSettings.biomartFolder 'biomartEnsemblGeneToGeneOntology.mat']);
        arrayDatabaseEnsemblGeneIDs = structBioMartData.arrayDatabaseEnsemblGeneIDs;
        arrayDatabaseGOTerms = structBioMartData.arrayDatabaseGOTerms;
        arrayDatabaseGONums = structBioMartData.arrayDatabaseGONums;
        arrayGOOfIntFlag = structBioMartData.arrayGOOfIntFlag;
        
    end
    
    
    numGenesMatchingGOTerm = sum(arrayGOOfIntFlag);
    
    arrayOutputRows = find(arrayGOOfIntFlag);
    
    arrayOutputEnsGeneIDs = zeros(numGenesMatchingGOTerm,1,'uint32');
    arrayOutputGOTerms = cell(numGenesMatchingGOTerm,1);
    arrayOutputGONums = zeros(numGenesMatchingGOTerm,1, 'uint32');
    
    for iOutputRow = 1:numGenesMatchingGOTerm,
        
        numDatabaseRow = arrayOutputRows(iOutputRow);
        
        arrayOutputEnsGeneIDs(iOutputRow) = arrayDatabaseEnsemblGeneIDs(numDatabaseRow);
        arrayOutputGOTerms{iOutputRow} = arrayDatabaseGOTerms{numDatabaseRow};
        arrayOutputGONums(iOutputRow) = arrayDatabaseGONums(numDatabaseRow);
        
    end
    
end

