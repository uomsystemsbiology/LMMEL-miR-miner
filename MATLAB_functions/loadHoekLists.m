function [ structCompLists ] = loadHoekLists( structSettings )
%LOADHOEKLISTS Summary of this function goes here
%   Detailed explanation goes here


    %%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
    % Load Various Gene Lists for Comparison with the Results
     %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 

    structCompLists = struct('type', cell(1), 'groupNames', cell(1), 'groupMembers', cell(1), 'stringLengths', cell(1));
 
    fileHoekList=fopen([ structSettings.InputFolder 'hoek_gene_list.csv' ],'rt');
    arrayTempHoekLists = textscan(fileHoekList,'%s','delimiter','\n','bufsize',4000000);
    fclose(fileHoekList);

    structCompLists.type = 'HoekGenes';
    %bring in group names from the header row
    stringHeaderRow = arrayTempHoekLists{1}{1};
    arrayCommaPos = strfind(stringHeaderRow, ',');
    structCompLists.groupNames = cell(2,1);
    structCompLists.groupNames{1} = stringHeaderRow(1:(arrayCommaPos(1)-1));
    structCompLists.groupNames{2} = stringHeaderRow((arrayCommaPos(1)+1):end); 
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
    structCompLists.stringLengths = cell(2,1);
    structCompLists.stringLengths{1} = zeros(numGrpOneMembs,1, 'int16');
    structCompLists.stringLengths{2} = zeros(numGrpTwoMembs,1, 'int16');
    structCompLists.groupMembers = cell(2,1);
    structCompLists.groupMembers{1} = cell(numGrpOneMembs,1);
    structCompLists.groupMembers{2} = cell(numGrpTwoMembs,1);
    for iRow = 2:length(arrayTempHoekLists{1}),
        stringRow = arrayTempHoekLists{1}{iRow};
        arrayCommaPos = strfind(stringRow, ',');
        strGrpOne = stringRow(1:(arrayCommaPos(1)-1));
        strGrpTwo = stringRow((arrayCommaPos(1)+1):end);
        if any(isspace(strGrpTwo)),
            strGrpTwo = strGrpTwo(~isspace(strGrpTwo));
        end
        if ~isempty(strGrpOne),
            structCompLists.groupMembers{1}{iRow-1} = strGrpOne;
            structCompLists.stringLengths{1}(iRow-1) = length(strGrpOne);
        end
        if ~isempty(strGrpTwo),
           structCompLists.groupMembers{2}{iRow-1} = strGrpTwo;
           structCompLists.stringLengths{2}(iRow-1) = length(strGrpTwo);
        end
    end
    
end

