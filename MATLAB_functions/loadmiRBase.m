function [ arrayTempMicRNAMatProdNames, arrayTempMicRNAMatProdAccNums, arrayTempMicRNAMatProdAltNames, arrayTempMicRNANames, arrayTempMicRNAAccNums ] = loadmiRBase( structFuncSettings )
%loadmiRBase :: A function to 
%   

%
% Created by Joe Cursons, 25/05/2015
% joseph (dot) cursons (at) unimelb (dot) edu (dot) au

%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 
% specify any fixed variables/folder paths
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 
    if isfield(structFuncSettings, 'miRBaseFolder'),
        stringDataFolder = structFuncSettings.miRBaseFolder;
    else
        stringDataFolder = 'C:\db\mirbase_21\';
    end
    
    if ~isfield(structFuncSettings, 'processMirBase'),
        %by default, don't overwrite previously extracted data
        structFuncSettings.processMirBase = false;
    end
    

    if structFuncSettings.processMirBase,
    
    %%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 
    % load the file in to memory
    %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %             
        %load the miRBase .dat file
        disp('Extracting micro-RNA data from the miRBase database');
        filemiRBase = fopen([stringDataFolder 'miRNA.dat']);
        arrayCellmiRBase = textscan(filemiRBase,'%s','delimiter','\n');
        arrayCellmiRBase = arrayCellmiRBase{1};
        fclose(filemiRBase);
        
        %determine the total number of entries and do some pre-indexing
        numTotalRows = length(arrayCellmiRBase);
        arrayHuEntryStartFlag = strncmp('ID   hsa', arrayCellmiRBase, 8);
        arrayHuEntryStartIndices = find(arrayHuEntryStartFlag);
        numHuMirBaseEntries = length(arrayHuEntryStartIndices);
        arrayEndEntryFlag = strncmp('//', arrayCellmiRBase, 2);
        arrayEndEntryIndices = find(arrayEndEntryFlag);
        
        arrayAccNumFlag = strncmp('AC', arrayCellmiRBase, 2);
        arrayAccNumIndices = find(arrayAccNumFlag);
        
        arrayFTFlag = strncmp('FT', arrayCellmiRBase, 2);
        arrayFTIndices = find(arrayFTFlag);

        %load the miRBase alias file
        filemiRBaseAliases = fopen([stringDataFolder 'aliases.txt']);
        arrayCellmiRBaseAliases = textscan(filemiRBaseAliases,'%s','delimiter','\n');
        arrayCellmiRBaseAliases = arrayCellmiRBaseAliases{1};
        fclose(filemiRBaseAliases);
        
        %miRBase does some painful sub-indexing system that works along the
        % lines of:
        
        % --> the script:
        %       1) pre-processes the file and identifies full entries which
        %           begin with "hsa-"
        %       2) looks at the "FT" lines of these entries, and pulls out
        %           the -3p and -5p products as separate entries
        %       3) moves through the alias file to extract information
        %           including "star" mapping


         %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
        %% define the output variables
         %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
        %assume that there are two mature products per entry
        
        arrayTempMicRNAMatProdNames = cell(numHuMirBaseEntries*2, 1);
        arrayTempMicRNAMatProdAccNums = cell(numHuMirBaseEntries*2, 1);
        arrayTempMicRNAMatProdAltNames = cell(numHuMirBaseEntries*2, 1);
        arrayTempMicRNANames = cell(numHuMirBaseEntries*2,1);
        arrayTempMicRNAAccNums = cell(numHuMirBaseEntries*2,1);
        arrayHasEntry = false(numHuMirBaseEntries*2,1);

        for iHuEntry = 1:numHuMirBaseEntries,
            
            numOutBaseIndex = (iHuEntry-1)*2;
        
            numStartRow = arrayHuEntryStartIndices(iHuEntry);
            numEndRowIndex = find(arrayEndEntryIndices > numStartRow, 1, 'first');
            numEndRow = arrayEndEntryIndices(numEndRowIndex);

            
            stringFirstRow = arrayCellmiRBase{numStartRow};
            arraySpacePos = strfind(stringFirstRow, ' ');
            
            %extract the miR name
            numNamePos = strfind(stringFirstRow, 'hsa-');
            numEndPosIndex = find(arraySpacePos > numNamePos, 1, 'first');
            numNameEndPos = arraySpacePos(numEndPosIndex) - 1;
            stringName = stringFirstRow(numNamePos:numNameEndPos);
                
            
            %extract the accession number
            arrayAccNumIndex = find(arrayAccNumIndices > numStartRow, 1, 'first');
            numAccNumRow = arrayAccNumIndices(arrayAccNumIndex);
            
            stringAccNumRow = arrayCellmiRBase{numAccNumRow};
            numAccNumPos = strfind(stringAccNumRow, 'MI');
            numSemiColonPos = strfind(stringAccNumRow, ';');
            
            numAccNums = length(numSemiColonPos);
            
            
            %find the FT rows for this entry
            arrayFTRowsForEntryIndices = (arrayFTIndices > numStartRow) & (arrayFTIndices < numEndRow);
            arrayFTRowsForEntryIndices = arrayFTIndices(arrayFTRowsForEntryIndices);
            numFTRows = length(arrayFTRowsForEntryIndices);
            arrayHasProductFlag = false(numFTRows,1);
            arrayHasAccessionFlag = false(numFTRows,1);
            for iRow = 1:numFTRows,
                numRow = arrayFTRowsForEntryIndices(iRow);
                stringRow = arrayCellmiRBase{numRow};
                if ~isempty(strfind(stringRow, '/product=')),
                    arrayHasProductFlag(iRow) = true;
                end
                if ~isempty(strfind(stringRow, '/accession=')),
                    arrayHasAccessionFlag(iRow) = true;
                end
            end
            numMatProds = sum(arrayHasProductFlag);
            arrayProdRowIndices = find(arrayHasProductFlag);
            arrayAccRowIndices = find(arrayHasAccessionFlag);
            for iProd = 1:numMatProds,
                
                stringProductRow = arrayCellmiRBase{arrayFTRowsForEntryIndices(arrayProdRowIndices(iProd))};
                stringMatProdAccessionRow = arrayCellmiRBase{arrayFTRowsForEntryIndices(arrayAccRowIndices(iProd))};
                
                arrayProdQuotationPos = strfind(stringProductRow, '"');
                arrayAccQuotationPos = strfind(stringMatProdAccessionRow, '"');
                
                stringMatProdAccNum = stringMatProdAccessionRow((arrayAccQuotationPos(1)+1):(arrayAccQuotationPos(2)-1));
                
                arrayTempMicRNAMatProdNames{numOutBaseIndex + iProd} = stringProductRow((arrayProdQuotationPos(1)+1):(arrayProdQuotationPos(2)-1));
                arrayTempMicRNAMatProdAccNums{numOutBaseIndex + iProd} = stringMatProdAccNum;
                arrayHasEntry(numOutBaseIndex + iProd) = true;
                

                %look for the accession number within the alias table
                arrayAliasMatchFlag = strncmp(stringMatProdAccNum, arrayCellmiRBaseAliases, length(stringMatProdAccNum));
                arrayAliasMatchIndex = find(arrayAliasMatchFlag);
                if ~isempty(arrayAliasMatchIndex),
                    stringAliasesRow = arrayCellmiRBaseAliases{arrayAliasMatchIndex};
                    arrayTabPos = strfind(stringAliasesRow, char(9));
                    arraySemiColonPos = strfind(stringAliasesRow, ';');
                    numAliases = length(arraySemiColonPos);

                    arrayTempMicRNAMatProdAltNames{numOutBaseIndex + iProd} = cell(numAliases, 1);
                    for iAlias = 1:numAliases,
                        if iAlias == 1,
                            numLastSpaceIndex = max(arrayTabPos);
                            numStartIndex = numLastSpaceIndex+1;
                        else
                            numLastSemiColonIndex = arraySemiColonPos(iAlias-1);
                            numStartIndex = numLastSemiColonIndex+1;
                        end

                        numEndIndex = arraySemiColonPos(iAlias)-1;

                        arrayTempMicRNAMatProdAltNames{numOutBaseIndex + iProd}{iAlias} = stringAliasesRow(numStartIndex:numEndIndex);

                    end

                end

                               
                
                arrayTempMicRNANames{numOutBaseIndex + iProd} = stringName;
                
                
                if numAccNums > 1,
                    disp(['warning: ' stringName ' has multiple miRBase accession numbers']);
                else
                    arrayTempMicRNAAccNums{numOutBaseIndex + iProd} = stringAccNumRow(numAccNumPos:(numSemiColonPos-1));
                end

            end
            
                 
        end
        
        numFinalEntries = sum(arrayHasEntry);
        arrayFinalEntryIndices = find(arrayHasEntry);
        
        
        arrayMicRNAMatProdNames = cell(numFinalEntries, 1);
        arrayMicRNAMatProdAccNums = cell(numFinalEntries, 1);
        arrayMicRNAMatProdAltNames = cell(numFinalEntries, 1);
        arrayMicRNANames = cell(numFinalEntries,1);
        arrayMicRNAAccNums = cell(numFinalEntries,1);
        for iOut = 1:numFinalEntries,
            arrayMicRNAMatProdNames{iOut} = arrayTempMicRNAMatProdNames{arrayFinalEntryIndices(iOut)};
            arrayMicRNAMatProdAccNums{iOut} = arrayTempMicRNAMatProdAccNums{arrayFinalEntryIndices(iOut)};
            arrayMicRNANames{iOut} = arrayTempMicRNANames{arrayFinalEntryIndices(iOut)};
            arrayMicRNAAccNums{iOut} = arrayTempMicRNAAccNums{arrayFinalEntryIndices(iOut)};
            
            arrayMicRNAMatProdAltNames{iOut} = cell(length(arrayTempMicRNAMatProdAltNames{arrayFinalEntryIndices(iOut)}),1);
            for iAltName = 1:length(arrayTempMicRNAMatProdAltNames{arrayFinalEntryIndices(iOut)}),
                arrayMicRNAMatProdAltNames{iOut}{iAltName} = arrayTempMicRNAMatProdAltNames{arrayFinalEntryIndices(iOut)}{iAltName};
            end
        end
        
        save([structFuncSettings.miRBaseFolder 'procMirBaseData.mat'], 'arrayMicRNAMatProdNames', 'arrayMicRNAMatProdAccNums', 'arrayMicRNAMatProdAltNames', 'arrayMicRNANames', 'arrayMicRNAAccNums');
        
    else
        
        structLoaded = load([structFuncSettings.miRBaseFolder 'procMirBaseData.mat']);
        arrayTempMicRNAMatProdNames = structLoaded.arrayMicRNAMatProdNames;
        arrayTempMicRNAMatProdAccNums = structLoaded.arrayMicRNAMatProdAccNums;
        arrayTempMicRNAMatProdAltNames = structLoaded.arrayMicRNAMatProdAltNames;
        arrayTempMicRNANames = structLoaded.arrayMicRNANames;
        arrayTempMicRNAAccNums = structLoaded.arrayMicRNAAccNums;
        
    end


end

