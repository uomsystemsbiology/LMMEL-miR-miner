function [ arrayUniqueGONums, arrayUniqueGenes, arrayGOMembershipMatrix ] = extractFullHumanGeneOntology( structFuncSettings )
%EXTRACTENTREZGENEINFO Summary of this function goes here
%   Detailed explanation goes here


    if structFuncSettings.performNewExtraction,

        disp([ char(9) '... extracting gene ontology information from the local database']);
        
        %open the input file 
        fileInputData = fopen([ structFuncSettings.dataFolder 'gene_association.goa_human' ],'rt');
        arrayInputCell = textscan(fileInputData,'%s','delimiter','\n');
        fclose(fileInputData);
        arrayGeneOntologyFile = arrayInputCell{1};
        clear arrayInputCell;

        numRows = length(arrayGeneOntologyFile);

        %move through all rows and extract all GO terms of interest and
        % corresponding gene names
        arrayAllGONums = zeros(numRows,1, 'uint32');
        arrayAllGenes = cell(numRows,1);
        for iRow = 1:numRows,
            stringRow = arrayGeneOntologyFile{iRow};
            %check for row exclusion marker
            if ~strncmp(stringRow, '!', 1),
                arrayTabPos = strfind(stringRow, char(9));
                %GO term falls within the fourth column
                stringGOTerm = stringRow((arrayTabPos(4)+1):(arrayTabPos(5)-1));
                %extract the numerical value for the GO term
                if strncmp(stringGOTerm, 'GO:', 3),
                    stringGONum = stringGOTerm(4:end);
                    numGOTerm = uint32(str2double(stringGONum));
                    arrayAllGONums(iRow) = numGOTerm;
                end
                                
                %HGNC falls within the third column
                stringHGNCSymbol = stringRow((arrayTabPos(2)+1):(arrayTabPos(3)-1));
                
                arrayAllGenes{iRow} = stringHGNCSymbol;
                
            else
                
                arrayAllGenes{iRow} = '';
            end
        end
        
        %create a unique list for both
        arrayUniqueGONums = unique(arrayAllGONums);
        arrayUniqueGenes = unique(arrayAllGenes);
        %extract the string length of the gene names for subsequent 
        % indexing
        arrayUniqueGeneStringLength = zeros(length(arrayUniqueGenes),1,'uint8');
        for iGene = 1:length(arrayUniqueGenes),
            arrayUniqueGeneStringLength(iGene) = length(arrayUniqueGenes{iGene});
        end
        
        %create the output membership matrix
        arrayGOMembershipMatrix = false(length(arrayUniqueGenes), length(arrayUniqueGONums));
        
        for iGOTerm = 1:length(arrayUniqueGONums),
            numGOTerm = arrayUniqueGONums(iGOTerm);
            
            arrayGOTermMatchIndices = arrayAllGONums == numGOTerm;
            
            arrayGenesWithThisGOTerm = arrayAllGenes(arrayGOTermMatchIndices);
            for iGene = 1:length(arrayGenesWithThisGOTerm),
                stringGene = arrayGenesWithThisGOTerm{iGene};
                
                numUniqueGeneIndex =  strncmp(stringGene,arrayUniqueGenes,length(stringGene)) & ...
                                           (arrayUniqueGeneStringLength == length(stringGene)) ;

                arrayGOMembershipMatrix(numUniqueGeneIndex,iGOTerm) = true;
                
            end
            
        end
        
        
        save([ structFuncSettings.dataFolder 'procFullHumanGeneOntology.mat'], 'arrayUniqueGONums', 'arrayUniqueGenes', 'arrayGOMembershipMatrix');
        
    else
        
        disp([ char(9) '... loading gene ontology information from a processed version of the local database']);
        
        structGeneOntoData = load([ structFuncSettings.dataFolder 'procFullHumanGeneOntology.mat']);
        arrayUniqueGONums = structGeneOntoData.arrayUniqueGONums;
        arrayUniqueGenes = structGeneOntoData.arrayUniqueGenes;
        arrayGOMembershipMatrix = structGeneOntoData.arrayGOMembershipMatrix;
                
    end
    
end

