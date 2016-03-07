function [ arrayOutputEnsGeneIDs, arrayOutputGOTerms, arrayOutputGONums ] = biomartEnsemblGeneToGeneOntology( arrayInputGOTerms, structFuncSettings )
%% [ arrayOutputEnsGeneIDs, arrayOutputGOTerms, arrayOutputGONums ] = biomartEnsemblGeneToGeneOntology( arrayInputGOTerms, structFuncSettings )
% This function is designed to take in a cell array of strings and match
%  them to specific GO categories. The function also has an input 
%  structured array specifying the location of the TargetScan data files,  
%  and a boolean flag to control re-processing of the original data.
%
%  Inputs:
%   - arrayInputGOTerms: a cell array of strings for searching the GO
%           categories.
%   - structSettings: a structured array which contains a number of
%           fields that control execution of this analysis. Fields required 
%           by this function include:
%       'biomartFolder' - a string specifying the folder path for the 
%           biomart
%       'performNewExtraction' - a Boolean flag specifying whether the
%           original data should be re-processed (True), or whether 
%           pre-processed data (saved as .mat) should be loaded (False). 
%           The script defaults to False as it is expected that input
%           strings will change between successive executions of the
%           function.
%   
%  Outputs:
%       A collection of cell arrays of strings, matched by row:
%           - arrayOutputEnsGeneIDs: miRBase name for the processed micro-RNA
%           - arrayOutputGOTerms: HGNC symbol for the mRNA
%           - arrayOutputGONums: RefSeq identifier for the mRNA
%       a uint32 vector, matched by row:
%           - arrayOutputEnsGeneIDs
%               
%  MATLAB Toolbox Dependencies:
%   - none
%
%  Function dependencies
%   - none
% 
% Users are encouraged to examine the Ensembl website and ensure that
%  the latest version of BioMart (currently build 83) is being used 
%  appropriately:
%       http://ensembl.org/biomart/
%
% For further information on Ensembl, please refer to:
%  Ensembl 2015
%  Fiona Cunningham, M. Ridwan Amode, Daniel Barrell, Kathryn Beal, et al.
%   [.. Stephen J. Trevanion, Andy Yates, Daniel R. Zerbino & Paul Flicek]
%  Nucleic Acids Research 2015 43 Database issue:D662-D669
%       http://dx.doi.org/10.1093/nar/gku1010
%
% For further information on Ensembl BioMart, please refer to: 
%  Ensembl BioMarts: a hub for data retrieval across taxonomic space.
%  Kinsella RJ, Kähäri A, Haider S, Zamora J, Proctor G, Spudich G,
%    Almeida-King J, Staines D, Derwent P, Kerhornou A, Kersey P, Flicek P. 
%  Database (Oxford). Vol. 2011 Published online Jul 23, 2011 
%       http://dx.doi.org/10.1093/database/bar030
%       PMID: 21785142
%
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 08/03/16
% 
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform variables/folder path pre-processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 
    %check the extraction/pre-processing setting has been assigned
    if ~isfield(structFuncSettings, 'performNewExtraction'),
        %by default just load the data, don't re-run the extraction
        structSettings.processTargetScan = false;
    end

    %check the dB folder path has been assigned
    if ~isfield(structFuncSettings, 'TargetScanFolder'),
        %default to the location on JCs desktop
        structSettings.TargetScanFolder = 'C:\db\targetscan_7p0\';
    end

         

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

