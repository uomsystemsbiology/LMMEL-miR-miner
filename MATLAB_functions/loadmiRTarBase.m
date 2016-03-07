function [ arrayMicroRNAIDs, arrayMicroRNANames, arrayMessRNANames, arrayMessRNAEntrezIDs, arraySpecies, arrayExperiment ] = loadmiRTarBase( structSettings )
%% [ arrayMicroRNAIDs, arrayMicroRNANames, arrayMessRNANames, arrayMessRNAEntrezIDs, arraySpecies, arrayExperiment ] = loadmiRTarBase( structSettings )
% This function is designed to take in a structured array which specifies
%  the location of the miRTarBase data files, and a boolean flag to control
%  re-processing of the original data.
%
%  Inputs:
%   - structSettings: a structured array which contains a number of
%           fields that control execution of different functions used
%           throughout this analysis. Fields required by this function
%           include:
%       'mirTarBaseFolder' - a string specifying the folder path for the 
%           miRTarBase data files
%       'processMirTarBase' - a Boolean flag specifying whether the
%           original data should be re-processed (True; doesn't take too 
%           long in this case as miRTarBase isn't *too* big), or whether 
%           pre-processed data (saved as .mat) is loaded (False; faster)
%
%  Output:
%       A collection of cell arrays of strings, matched by row:
%           - arrayMicroRNANames: miRBase name for the processed micro-RNA
%           - arrayMessRNANames: HGNC symbol for the mRNA
%           - arraySpecies: Species for the miR and mRNA (separated by
%               '//')
%           - arrayExperiment: List ('//' separator) of experimental
%               evidence which supports the relationship
%       and uint32 vectors, matched by row:
%           - arrayMicroRNAIDs: miRBase identifier for the micro-RNA (MIRT
%               prefix)
%           - arrayMessRNAEntrezIDs: Entrez identifier for the mRNA
%               
%  MATLAB Toolbox Dependencies:
%   - none
%
%  Function dependencies
%   - none
%
% Users are encouraged to examine the miRTarBase website and ensure that
%  the latest version (currently 6.1) is being used appropriately:
%   http://mirtarbase.mbc.nctu.edu.tw/
% For further information on miRTarBase, please refer to:
% miRTarBase 2016: updates to the experimentally validated miRNA-target
%   interactions database. 
%  Chou CH, Chang NW, Shrestha S, Hsu SD, Lin YL, Lee WH et al. 
%  Nucleic acids research. 2016; 44(D1): D239-47. 
%  http://dx.doi.org/10.1093/nar/gkv1258
% 
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 04/03/16
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform pre-processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %check the dB folder path has been assigned
    if ~isfield(structSettings, 'mirTarBaseFolder'),
        %default to the file path on JCs work PC
        structSettings.mirTarBaseFolder = 'C:\db\miRTarBase_6p1\';
    end
    %check the extraction/pre-processing setting has been assigned
    if ~isfield(structSettings, 'processMirTarBase'),
        %default to re-processing the original data
        structSettings.processMirTarBase = true;
    end


    if structSettings.processMirTarBase,
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Process the data and populate the output arrays
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %      

        %determine the location of the miRTarBase database file
        stringFilePath = [ structSettings.mirTarBaseFolder 'hsa_MTI.csv' ];

        %load the data file into a cell array
        filemiRTarBase = fopen(stringFilePath);
        arrayFile = textscan( filemiRTarBase, '%s', 'delimiter', '\n' );
        fclose(filemiRTarBase);
        arraymiRTarBase = arrayFile{1};

        %examine the size of the array
        numFileLines = length(arraymiRTarBase);
        
        %ignore the header reow
        arrayOutputLines = 2:numFileLines;
        numOutputLines = numFileLines - 1;

        %create output arrays
        arrayMicroRNAIDs = zeros(numOutputLines,1,'uint32');
        arrayMicroRNANames = cell(numOutputLines,1);
        arrayMessRNANames = cell(numOutputLines,1);
        arrayMessRNAEntrezIDs = zeros(numOutputLines,1,'uint32');
        arraySpecies = cell(numOutputLines,1);
        arrayExperiment = cell(numOutputLines,1);
         
        
        %progress through the miRTarBase data file
        for iLine = arrayOutputLines,
            
            %determine the corresponding output row
            numOutRow = iLine - 1;

            %extract the line for this row
            stringLine = arraymiRTarBase{iLine};
            %determine the location of separators (,)
            arrayCommaPos = strfind(stringLine, ',');
            
            %extract the miR ID from the first column
            stringMicRNAID = stringLine(1:(arrayCommaPos(1)-1));
            %check the prefix is as expected
            if strncmp(stringMicRNAID(1:4), 'MIRT', 4),
                %convert to double
                numMicRNAIDDouble = str2double(stringMicRNAID(5:end));
                %and check that it is a valid uint32 value for output
                if (numMicRNAIDDouble < ((2^32) - 1)),
                    arrayMicroRNAIDs(numOutRow) = uint32(numMicRNAIDDouble);
                else
                    disp('warning: miRTarBase IDs are too large to store as uint32');
                end
            else
                disp('warning: mirTarBase ID does not have a MIRT prefix');
            end
            
            %extract the miR name in the second column
            stringMicRNAName = stringLine((arrayCommaPos(1)+1):(arrayCommaPos(2)-1));
            arrayMicroRNANames{numOutRow} = stringMicRNAName;
            
            %extract the species in the third column and sixth column
            stringSpeciesOne = stringLine((arrayCommaPos(2)+1):(arrayCommaPos(3)-1));
            stringSpeciesTwo = stringLine((arrayCommaPos(5)+1):(arrayCommaPos(6)-1));
            if strncmp(stringSpeciesOne, stringSpeciesTwo, length(stringSpeciesTwo)),
                arraySpecies{numOutRow} = stringSpeciesOne;
            else
                arraySpecies{numOutRow} = [ stringSpeciesOne '//' stringSpeciesTwo ];
            end
            
            
            %extract the target name in the fourth column
            stringMessRNAName = stringLine((arrayCommaPos(3)+1):(arrayCommaPos(4)-1));
            arrayMessRNANames{numOutRow} = stringMessRNAName;
            
            %extract the target ID in the fifth column
            stringMessRNAID = stringLine((arrayCommaPos(4)+1):(arrayCommaPos(5)-1));
            %convert to a double
            numMessRNAIDDouble = str2double(stringMessRNAID);
            %check whether the number is valid for a uint32 integer
            if (numMessRNAIDDouble < ((2^32) - 1)),
                arrayMessRNAEntrezIDs(numOutRow) = uint32(numMessRNAIDDouble);
            else
                if ~isnan(numMessRNAIDDouble),
                    disp('warning: mRNA Entrez ID is too large to store as uint32');
                end
            end
            
            %extract the experiment in the seventh column
            stringExperiments = stringLine((arrayCommaPos(6)+1):(arrayCommaPos(7)-1));
            arrayExperiment{numOutRow} = stringExperiments;
            
        end
         
        %save the output vectors to a .mat file so that they can be loaded
        % in the future
        save([ structSettings.mirTarBaseFolder 'PreProcmiRTarBase.mat' ], 'arrayMicroRNAIDs', 'arrayMicroRNANames', 'arrayMessRNANames', 'arrayMessRNAEntrezIDs', 'arraySpecies', 'arrayExperiment');

    else
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Load the .mat file and populate the output vectors
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %     
        
        %load the .mat file, which stores each vector within a structured
        % array
        structLoaded = load([ structSettings.mirTarBaseFolder 'PreProcmiRTarBase.mat' ]);
        
        %assign the data to the correct output vector
        arrayMicroRNAIDs = structLoaded.arrayMicroRNAIDs;
        arrayMicroRNANames = structLoaded.arrayMicroRNANames;
        arrayMessRNANames = structLoaded.arrayMessRNANames;
        arrayMessRNAEntrezIDs = structLoaded.arrayMessRNAEntrezIDs;
        arraySpecies = structLoaded.arraySpecies;
        arrayExperiment = structLoaded.arrayExperiment;
        
    end
    
    
end

