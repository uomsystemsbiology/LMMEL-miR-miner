function [ structOutputData, structClusters ] = loadLudwigData( structInputSettings )
%% [ structOutputData, structClusters ] = loadLudwigData( structSettings )
% This function is designed to take in a structured array which specifies
%  the location of the LM-MEL panel data and extract the miR and mRNA data 
%  into structured arrays with fields for meta-information extracted from
%  the pre-processed files (HGNC symbol; RefSeq identifier etc).
%
%  Inputs:
%   - structInputSettings: a structured array which contains a number of
%           fields that control execution of different functions used
%           throughout this analysis. Fields required by this function
%           include:
%       'InputFolder' - a string specifying the folder path for the LM-MEL
%           panel processed data files
%
%  Output:
%   - structOutputData: a 1D, length 2, structured array which contains
%           fields that specify information for:
%               (1) miR data
%               (2) mRNA data
%           These fields include:
%               'DataType' - string specifying whether the structured array 
%                   contains miR or mRNA data
%               'CellLine' -  cell array of strings containing the LM-MEL 
%                   cell line number/labelling
%               'Target' - cell array of strings containing the HGNC (mRNA)
%                   or miRBase (miR) label for the RNA transcript
%               'TargetRefSeq' - cell array of strings containing the 
%                   RefSeq accession number for mRNA transcripts
%               'TargetID' - 
%               'TargetEntrez' - 
%               'TargetILMN' - 
%               'Data' - 2D double-precision array with RNA transcript abundance,
%                   indexed by cell line and target
%   - structClusters: a 1D, length 2, structured array which contains 
%           descriptive fields for clustering by:
%               (1) pigmentation status (NB: clustering by mRNA expression)
%               (2) invasiveness status (independent matrigel assay)
%           The fields include:
%               'type' - a string specifying the type of data used for
%                   clustering
%               'groupNames' - cell array of strings specifying the
%                   descriptor (high/low) associated with each
%                   group/cluster
%               'groupMembers' - cell array of strings speciying the cell
%                   lines assopciated with each group/cluster
%
%  MATLAB Toolbox Dependencies:
%   - ?
%
% This MATLAB function has been released for the manuscript which is under 
%   review at BMC Genome Biology:
%   M.C. Andrews, J. Cursons et al. (2016). Systems analysis of the Ludwig 
%       Melbourne melanoma cell line panel identifies miR-29b repression as
%        a driver of melanoma invasiveness.
%   doi: not-yet-assigned
% 
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 03/03/16
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% User-specified settings
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %specify the filenames for the processed miR and mRNA data
    arrayFileNames = { [ 'LM-MEL-57 miRNA Transpose.csv' ];
                       [ 'LM-MEL-57 GE Transpose Annot.csv' ] };

    %specify the expected (known) size of the data because it's a lot 
    % quicker when creating the output arrays
    arrayExpectedSize = { [ 2593 58 ];      %rows cols
                          [ 47232 62 ] };

    %specify the number of header columns in each processed data file
    arrayNumHeaderCols = { [1];
                           [5] };                   
    %and the number of header rows  
    arrayNumHeaderRows = { [1];
                           [1] };
      
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Load/Extract the data
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %               

    %check the input settings
    if isempty(structInputSettings),
        disp('warning: the input settings for the ludwig data structured array has not been properly defined');
    else
        %check the the data folder path has been specified
        if isempty(structInputSettings.InputFolder),
            %default to the known path (note that this will not work on
            % other machines without being modified)
            structInputSettings.InputFolder = 'C:\wc\2015_ludwig_melanoma\data\';
        end
    end
    
    %read the processed data from each file into large cell arrays
    arrayContent = cell(length(arrayFileNames),1);
    for iFile = 1:length(arrayFileNames),
        fid=fopen([ structInputSettings.InputFolder arrayFileNames{iFile} ],'rt');
        arrayContent{iFile} = textscan(fid,'%s','delimiter','\n');
        fclose(fid);
    end
    
    %read in the pigmentation clustering information
    filePigmentClusters = fopen([structInputSettings.InputFolder 'ClusterCellLinesByPigment.txt']);
    arrayPigmentClusters = textscan(filePigmentClusters,'%s','delimiter','\n');
    fclose(filePigmentClusters);

    %read in the invasiveness clustering information
    fileInvasiveClusters = fopen([structInputSettings.InputFolder 'matrigel_invasiveness.csv']);
    arrayInvasiveClusters = textscan(fileInvasiveClusters,'%s','delimiter','\n');
    fclose(fileInvasiveClusters);
    
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform processing and populate output data arrays
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %           
 
    %create the output struct for the data
    structOutputData = struct( 'DataType', cell(length(arrayFileNames),1), ...
                               'CellLine', cell(length(arrayFileNames),1), ...
                               'Target', cell(length(arrayFileNames),1), ...
                               'TargetRefSeq', cell(length(arrayFileNames),1), ...
                               'TargetID', cell(length(arrayFileNames),1), ...
                               'TargetEntrez', cell(length(arrayFileNames),1), ...
                               'TargetILMN',cell(length(arrayFileNames),1), ...
                               'Data', cell(length(arrayFileNames),1));
    
    %move through each file
    for iFile = 1:length(arrayFileNames),
        
        %define the data type
        if strcmp(arrayFileNames{iFile}, 'LM-MEL-57 miRNA Transpose.csv'),
            structOutputData(iFile).DataType = 'miRNA';
        elseif strcmp(arrayFileNames{iFile}, 'LM-MEL-57 GE Transpose Annot.csv'),
            structOutputData(iFile).DataType = 'mRNA:MicroArray';
        end
        
        %extract the cell line data from the header row
        arrayHeaderRowString = arrayContent{iFile}{1}{1};
        arrayHeaderRowComma = strfind(arrayHeaderRowString,',');
        %check the input data integrity - # of columns
        if (length(arrayHeaderRowComma)+1) == arrayExpectedSize{iFile}(2),
            %if the data look OK, read in cell lines
            numCellLines = arrayExpectedSize{iFile}(2) - arrayNumHeaderCols{iFile};
            structOutputData(iFile).CellLine = cell(numCellLines,1);
            for iCellLine = 1:numCellLines,
                numStartIndex = arrayHeaderRowComma(iCellLine + (arrayNumHeaderCols{iFile}-1))+1;
                %tricksy csv
                if iCellLine < numCellLines,
                    numEndIndex = arrayHeaderRowComma(iCellLine+arrayNumHeaderCols{iFile})-1;
                else
                    numEndIndex = length(arrayHeaderRowString);
                end
                structOutputData(iFile).CellLine{iCellLine} = arrayHeaderRowString(numStartIndex:numEndIndex);
            end
        else
            disp('check that the data is reading in properly - unexepcted number of columns');
        end
        
        %check the input data integrity - # of rows
        if length(arrayContent{iFile}{1}) == arrayExpectedSize{iFile}(1),
            %if the data look OK, read in target info and data
            numTargets = length(arrayContent{iFile}{1}) - arrayNumHeaderRows{iFile};
            structOutputData(iFile).Target = cell(numTargets,1);
            structOutputData(iFile).Data = zeros(numTargets,numCellLines,'double');
            if ~isempty(strfind(structOutputData(iFile).DataType, 'mRNA'));
                structOutputData(iFile).TargetRefSeq = cell(numTargets,1);
                structOutputData(iFile).TargetID = cell(numTargets,1);
                structOutputData(iFile).TargetEntrez = cell(numTargets,1);
                structOutputData(iFile).TargetILMN = cell(numTargets,1);
            end
            
            for iRow = 2:length(arrayContent{iFile}{1}),
                arrayRowString = arrayContent{iFile}{1}{iRow};
                arrayCommaPos = strfind(arrayRowString,',');
                
                numTarget = iRow-1;
                if iFile == 1,
                    %reading the miRNA data, only one header column
                    structOutputData(iFile).Target{numTarget} = arrayRowString(1:(arrayCommaPos(1)-1));
                elseif iFile == 2,
                    %reading the miRNA data, several one header column
                    structOutputData(iFile).TargetID{numTarget} = arrayRowString(1:(arrayCommaPos(1)-1));
                    structOutputData(iFile).TargetEntrez{numTarget} = arrayRowString((arrayCommaPos(1)+1):(arrayCommaPos(2)-1));
                    structOutputData(iFile).TargetILMN{numTarget} = arrayRowString((arrayCommaPos(2)+1):(arrayCommaPos(3)-1));
                    structOutputData(iFile).TargetRefSeq{numTarget} = arrayRowString((arrayCommaPos(3)+1):(arrayCommaPos(4)-1));
                    structOutputData(iFile).Target{numTarget} = arrayRowString((arrayCommaPos(4)+1):(arrayCommaPos(5)-1));
                else
                    disp('The number of input rows does not match the expected value');
                end
                
                %read in all of the data
                for iCellLine = 1:numCellLines,
                    numCol = arrayNumHeaderCols{iFile}+iCellLine;
                    if iCellLine < numCellLines,
                        structOutputData(iFile).Data(numTarget,iCellLine) = str2double(arrayRowString((arrayCommaPos(numCol-1)+1):(arrayCommaPos(numCol)-1)));
                    else
                        structOutputData(iFile).Data(numTarget,iCellLine) = str2double(arrayRowString((arrayCommaPos(numCol-1)+1):end));
                    end
                end
                
            end

        else
            disp('check that the data is reading in properly - unexpected number of columns');
        end
        
    end
    
 
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform processing and populate output clustering arrays
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %        
    
    %create an output structure for the clustering information
    structClusters = struct('type', cell(2,1), 'groupNames', cell(2,1), 'groupMembers', cell(2,1)); 
 
    
    %clustering by pigmentation status
    structClusters(1).type = 'pigment';
    structClusters(1).groupNames = cell(2,1);
    structClusters(1).groupNames{1} = 'Low';
    structClusters(1).groupNames{2} = 'High';
    structClusters(1).groupMembers = cell(2,1);
    
    %count the number of cells associated with each pigment class to create the array
    numLowPigment = 0;
    numHighPigment = 0;
    for iRow = 2:length(arrayPigmentClusters{1}),
        if ~isempty(strfind(arrayPigmentClusters{1}{iRow}, 'Low')),
            numLowPigment = numLowPigment+1;
        elseif ~isempty(strfind(arrayPigmentClusters{1}{iRow}, 'High')),
            numHighPigment = numHighPigment+1;
        end
    end
    structClusters(1).groupMembers{1} = cell(numLowPigment,1);
    structClusters(1).groupMembers{2} = cell(numHighPigment,1);
    %populate the output array
    numLowPigment = 1;
    numHighPigment = 1;
    for iRow = 2:length(arrayPigmentClusters{1}),
        stringRow = arrayPigmentClusters{1}{iRow};
        arrayRowTabs = strfind(stringRow, char(9));
        stringCellLine = stringRow(1:(arrayRowTabs(1)-1));
        if ~isempty(strfind(arrayPigmentClusters{1}{iRow}, 'Low')),
            structClusters(1).groupMembers{1}{numLowPigment} = stringCellLine;
            numLowPigment = numLowPigment+1;
        elseif ~isempty(strfind(arrayPigmentClusters{1}{iRow}, 'High')),
            structClusters(1).groupMembers{2}{numHighPigment} = stringCellLine;
            numHighPigment = numHighPigment+1;
        end
    end


    
    %clustering by invasiveness status
    structClusters(2).type = 'invasiveness';
    structClusters(2).groupNames = cell(2,1);
    structClusters(2).groupNames{1} = 'Low';
    structClusters(2).groupNames{2} = 'High';
    structClusters(2).groupMembers = cell(2,1);
    
    %count the number of cells associated with each class to create the array
    numLowInvasive = 0;
    numHighInvasive = 0;
    for iRow = 2:length(arrayInvasiveClusters{1}),
        %process the row/string
        stringRow = arrayInvasiveClusters{1}{iRow};
        arrayCommaPos = strfind(stringRow, ',');
        stringHighInv = stringRow(1:(arrayCommaPos(1)-1));
        stringLowInv = stringRow((arrayCommaPos(1)+1):end);

        if ~isempty(stringLowInv),
            numLowInvasive = numLowInvasive+1;
        end
        if ~isempty(stringHighInv),
            numHighInvasive = numHighInvasive+1;
        end
    end
    structClusters(2).groupMembers{1} = cell(numLowInvasive,1);
    structClusters(2).groupMembers{2} = cell(numHighInvasive,1);
    %populate the output array
    numLowInvasive = 1;
    numHighInvasive = 1;
    for iRow = 2:length(arrayInvasiveClusters{1}),
        %process the row/string
        stringRow = arrayInvasiveClusters{1}{iRow};
        arrayCommaPos = strfind(stringRow, ',');
        stringHighInv = stringRow(1:(arrayCommaPos(1)-1));
        stringLowInv = stringRow((arrayCommaPos(1)+1):end);

        if ~isempty(stringLowInv),
            structClusters(2).groupMembers{1}{numLowInvasive} = stringLowInv;
            numLowInvasive = numLowInvasive+1;
        end
        if ~isempty(stringHighInv),
            structClusters(2).groupMembers{2}{numHighInvasive} = stringHighInv;
            numHighInvasive = numHighInvasive+1;
        end
    end
    
return

