function [ structOutputData, structClusters ] = loadLudwigData( structSettings )
%loadLudwigData :: Load the Ludwig melanoma data .csv for the
%   specified cell line into a structured array
%
% This function takes a specified cell line
% Inputs:
%   - structSettings: a structured array containing settings to control
%       output
% Outputs:
%   - structOutputData: a structured array containing miRNA and mRNA
%       abundance data across the different cell lines

arrayFileNames = { [ 'LM-MEL-57 miRNA Transpose.csv' ];
                   [ 'LM-MEL-57 GE Transpose Annot.csv' ] };
               
arrayExpectedSize = { [ 2593 58 ];  %row col
                      [ 47232 62 ] };

arrayNumHeaderCols = { [1];
                       [5] };                   
  
arrayNumHeaderRows = { [1];
                       [1] };
                   
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %  
   %%% Check Input
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %      
    if isempty(structSettings),
        % *** DO NOTHING ***
    else
        % *** Error ***
    end
    
    
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %  
   %%% Define paths
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %  
    stringCurrDir = cd;
        
    
    
    %create an output structure for clustering information
    structClusters = struct('type', cell(2,1), 'groupNames', cell(2,1), 'groupMembers', cell(2,1)); 

    
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %  
    %% Load the data array into MATLAB
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %       
    %read the data from each file into large cell arrays
    arrayContent = cell(length(arrayFileNames),1);
    for iFile = 1:length(arrayFileNames),
        fid=fopen([ structSettings.InputFolder arrayFileNames{iFile} ],'rt');
        arrayContent{iFile} = textscan(fid,'%s','delimiter','\n');
        fclose(fid);
    end
    
    
    %read in the pigmentation clustering information
    filePigmentClusters = fopen([structSettings.InputFolder 'ClusterCellLinesByPigment.txt']);
    arrayPigmentClusters = textscan(filePigmentClusters,'%s','delimiter','\n');
    fclose(filePigmentClusters);

    %read in the invasiveness clustering information
    fileInvasiveClusters = fopen([structSettings.InputFolder 'matrigel_invasiveness.csv']);
    arrayInvasiveClusters = textscan(fileInvasiveClusters,'%s','delimiter','\n');
    fclose(fileInvasiveClusters);
    
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %  
    %% Extract the required data
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   

    %create the output struct
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
    
       
       
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %  
    %% Load in Various Clustering Information
%   %%%%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   
    
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

