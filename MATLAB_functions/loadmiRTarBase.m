function [ arrayMicroRNAIDs, arrayMicroRNANames, arrayMessRNANames, arrayMessRNAEntrezIDs, arraySpecies, arrayExperiment ] = loadmiRTarBase( structSettings )
%LOADMIRBASE Summary of this function goes here
%   Detailed explanation goes here

     %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
    %% check input settings 
     %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
    %check the dB folder path
    if ~isfield(structSettings, 'mirTarBaseFolder'),
        %default to the file path on my work PC
        structSettings.mirTarBaseFolder = 'C:\db\miRTarBase_4p5\';
    end
    %and the extraction settings
    if ~isfield(structSettings, 'processMirTarBase'),
        %default to the file path on my work PC
        structSettings.processMirTarBase = true;
    end
    

    if structSettings.processMirTarBase,
         %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
        %% load the database file
         %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %

        %location of the miRTarBase database file
        stringFilePath = [ structSettings.mirTarBaseFolder 'hsa_MTI.csv' ];

        %load the TargetScan miR family file
        filemiRTarBase = fopen(stringFilePath);
        arrayFile = textscan( filemiRTarBase, '%s', 'delimiter', '\n' );
        fclose(filemiRTarBase);
        arraymiRTarBase = arrayFile{1};

        numFileLines = length(arraymiRTarBase);
        
        %we don't export the header
        arrayOutputLines = 2:numFileLines;
        numOutputLines = numFileLines - 1;

         %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
        %% create the output arrays
         %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
        arrayMicroRNAIDs = zeros(numOutputLines,1,'uint32');
        arrayMicroRNANames = cell(numOutputLines,1);
        arrayMessRNANames = cell(numOutputLines,1);
        arrayMessRNAEntrezIDs = zeros(numOutputLines,1,'uint32');
        arraySpecies = cell(numOutputLines,1);
        arrayExperiment = cell(numOutputLines,1);
         
        
         %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
        %% populate the output arrays
         %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
        for iLine = arrayOutputLines,
            
            numOutRow = iLine - 1;

            stringLine = arraymiRTarBase{iLine};
            arrayCommaPos = strfind(stringLine, ',');
            
            %miR ID in the first column
            stringMicRNAID = stringLine(1:(arrayCommaPos(1)-1));
            if strncmp(stringMicRNAID(1:4), 'MIRT', 4),
                numMicRNAIDDouble = str2double(stringMicRNAID(5:end));
                if (numMicRNAIDDouble < ((2^32) - 1)),
                    arrayMicroRNAIDs(numOutRow) = uint32(numMicRNAIDDouble);
                else
                    disp('warning: miRTarBase IDs are too large to store as uint32');
                end
            else
                disp('warning: mirTarBase ID does not have a MIRT prefix');
            end
            
            %miR name in the second column
            stringMicRNAName = stringLine((arrayCommaPos(1)+1):(arrayCommaPos(2)-1));
            arrayMicroRNANames{numOutRow} = stringMicRNAName;
            
            %species in the third column and sixth column
            stringSpeciesOne = stringLine((arrayCommaPos(2)+1):(arrayCommaPos(3)-1));
            stringSpeciesTwo = stringLine((arrayCommaPos(5)+1):(arrayCommaPos(6)-1));
            if strncmp(stringSpeciesOne, stringSpeciesTwo, length(stringSpeciesTwo)),
                arraySpecies{numOutRow} = stringSpeciesOne;
            else
                arraySpecies{numOutRow} = [ stringSpeciesOne '//' stringSpeciesTwo ];
            end
            
            
            %target name in the fourth column
            stringMessRNAName = stringLine((arrayCommaPos(3)+1):(arrayCommaPos(4)-1));
            arrayMessRNANames{numOutRow} = stringMessRNAName;
            
            %target ID in the fifth column
            stringMessRNAID = stringLine((arrayCommaPos(4)+1):(arrayCommaPos(5)-1));
            numMessRNAIDDouble = str2double(stringMessRNAID);
            if (numMessRNAIDDouble < ((2^32) - 1)),
                arrayMessRNAEntrezIDs(numOutRow) = uint32(numMessRNAIDDouble);
            else
                if ~isnan(numMessRNAIDDouble),
                    disp('warning: mRNA Entrez ID is too large to store as uint32');
                end
            end
            
            
            %experiment in the seventh column
            stringExperiments = stringLine((arrayCommaPos(6)+1):(arrayCommaPos(7)-1));
            arrayExperiment{numOutRow} = stringExperiments;
            
        end
         
        
        save([ structSettings.mirTarBaseFolder 'PreProcmiRTarBase.mat' ], 'arrayMicroRNAIDs', 'arrayMicroRNANames', 'arrayMessRNANames', 'arrayMessRNAEntrezIDs', 'arraySpecies', 'arrayExperiment');

    else
        
        arrayLoaded = load([ structSettings.mirTarBaseFolder 'PreProcmiRTarBase.mat' ]);
        arrayMicroRNAIDs = arrayLoaded.arrayMicroRNAIDs;
        arrayMicroRNANames = arrayLoaded.arrayMicroRNANames;
        arrayMessRNANames = arrayLoaded.arrayMessRNANames;
        arrayMessRNAEntrezIDs = arrayLoaded.arrayMessRNAEntrezIDs;
        arraySpecies = arrayLoaded.arraySpecies;
        arrayExperiment = arrayLoaded.arrayExperiment;
        
    end
    
    
end

