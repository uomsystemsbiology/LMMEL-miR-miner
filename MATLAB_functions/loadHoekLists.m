function [ structCompLists ] = loadHoekLists( structInputSettings )
%% [ structCompLists ] = loadHoekLists( structInputSettings )
% This function is designed to take in a structured array which specifies
%  the location of the Hoek proliferative/invasive gene lists and extract  
%  these into structured arrays.
%
%  Inputs:
%   - structInputSettings: a structured array which contains a number of
%           fields that control execution of different functions used
%           throughout this analysis. Fields required by this function
%           include:
%       'InputFolder' - a string specifying the folder path for the Hoek
%           proliferative/invasive gene lists
%
%  Output:
%   - structCompLists: a 1D, length 1, structured array which contains
%           fields that specify information for:
%               (1) signatures from Hoek (2006) - details below
%               (2) signatures from Widmer (2012) - details below
%           These fields include:
%               'type'
%               'groupNames'
%               'groupMembers'
%               'stringLengths'
%               
%  MATLAB Toolbox Dependencies:
%   - ?
%
% The gene lists extracted by this script were identified in the following
% manuscripts:
%
% KS Hoek et al. [D Schadendorf, R Dummer]. Metastatic potential of 
%       melanomas defined by specific gene expression profiles with no 
%       BRAF signature. Pigment Cell Res. 2006 Aug;19(4):290-302.
%           http://dx.doi.org/10.1111/j.1600-0749.2006.00322.x
%
% DS Widmer et al. [R Dummer, KS Hoek]. Systematic classification of 
%       melanoma cells by phenotype-specific gene expression mapping.
%       Pigment Cell Melanoma Res. 2012 May;25(3):343-53.
%           http://dx.doi.org/10.1111/j.1755-148X.2012.00986.x
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

    %specify the filenames containing the data
    arrayFileNames = { ['hoek_105genelist.csv'];
                       ['hoek_97genelist.csv'] };
    numFiles = length(arrayFileNames);

    %specify the strings for the structured array labels
    arrayListLabels = { ['Hoek 2006 (105 genes)'];
                        ['Widmer 2012 (97 genes)'] };
                   
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Load and Process the Data
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
                   
    %create the output structured array
    structCompLists = struct( 'type', cell(2,1), ...
                              'groupNames', cell(2,1), ...
                              'groupMembers', cell(2,1), ...
                              'stringLengths', cell(2,1) );
 
    %move through each of the specified files and extract the gene lists
    for iFile = 1:numFiles,
        
        %label the list within the structured array
        structCompLists(iFile).type = arrayListLabels{iFile};
        
        %load the gene list csv
        fileHoekList=fopen([ structInputSettings.InputFolder arrayFileNames{iFile} ],'rt');
        arrayTempHoekLists = textscan(fileHoekList,'%s','delimiter','\n');
        fclose(fileHoekList);

        %bring in group names from the header row
        stringHeaderRow = arrayTempHoekLists{1}{1};
        arrayCommaPos = strfind(stringHeaderRow, ',');
        structCompLists(iFile).groupNames = cell(2,1);
        structCompLists(iFile).groupNames{1} = stringHeaderRow(1:(arrayCommaPos(1)-1));
        structCompLists(iFile).groupNames{2} = stringHeaderRow((arrayCommaPos(1)+1):end); 
        %determine how many members are in each group
        numGrpOneMembs = 0;
        numGrpTwoMembs = 0;
        for iRow = 2:length(arrayTempHoekLists{1}),
            stringRow = arrayTempHoekLists{1}{iRow};
            arrayCommaPos = strfind(stringRow, ',');
            strGrpOne = stringRow(1:(arrayCommaPos(1)-1));
            strGrpTwo = stringRow((arrayCommaPos(1)+1):end);
            if any(isspace(strGrpTwo)),
                strGrpTwo = strGrpTwo(~isspace(strGrpTwo));
            end
            if ~isempty(strGrpOne),
                numGrpOneMembs = numGrpOneMembs+1;
            end
            if ~isempty(strGrpTwo),
                numGrpTwoMembs = numGrpTwoMembs+1;
            end
        end
        
        %populate the group members from the other rows 
        structCompLists(iFile).stringLengths = cell(2,1);
        structCompLists(iFile).stringLengths{1} = zeros(numGrpOneMembs,1, 'int16');
        structCompLists(iFile).stringLengths{2} = zeros(numGrpTwoMembs,1, 'int16');
        structCompLists(iFile).groupMembers = cell(2,1);
        structCompLists(iFile).groupMembers{1} = cell(numGrpOneMembs,1);
        structCompLists(iFile).groupMembers{2} = cell(numGrpTwoMembs,1);
        for iRow = 2:length(arrayTempHoekLists{1}),
            stringRow = arrayTempHoekLists{1}{iRow};
            arrayCommaPos = strfind(stringRow, ',');
            strGrpOne = stringRow(1:(arrayCommaPos(1)-1));
            strGrpTwo = stringRow((arrayCommaPos(1)+1):end);
            if any(isspace(strGrpTwo)),
                strGrpTwo = strGrpTwo(~isspace(strGrpTwo));
            end
            if ~isempty(strGrpOne),
                structCompLists(iFile).groupMembers{1}{iRow-1} = strGrpOne;
                structCompLists(iFile).stringLengths{1}(iRow-1) = length(strGrpOne);
            end
            if ~isempty(strGrpTwo),
               structCompLists(iFile).groupMembers{2}{iRow-1} = strGrpTwo;
               structCompLists(iFile).stringLengths{2}(iRow-1) = length(strGrpTwo);
            end
        end
    end
    
end

