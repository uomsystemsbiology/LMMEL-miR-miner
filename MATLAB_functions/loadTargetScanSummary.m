function [ arrayOutMicRNAName, arrayOutMessRNA, arrayOutMessRNARefSeq, arrayOutContPlus ] = loadTargetScanSummary( structSettings )
%% [ arrayOutMicRNAName, arrayOutMessRNA, arrayOutMessRNARefSeq, arrayOutContPlus ] = loadTargetScanSummary( structSettings )
% This function is designed to take in a structured array which specifies
%  the location of the TargetScan data files, and a boolean flag to 
%  control re-processing of the original data.
% NB: the context+ which is output from this function is rescaled to fit
%       into a uint32 format. In this case -1000 (0 -> -x rescaled to
%       0 -> x000).
%
%  Inputs:
%   - structSettings: a structured array which contains a number of
%           fields that control execution of different functions used
%           throughout this analysis. Fields required by this function
%           include:
%       'TargetScanFolder' - a string specifying the folder path for the 
%           TargetScan data files
%       'processTargetScan' - a Boolean flag specifying whether the
%           original data should be re-processed (True; takes a long time),
%           or whether pre-processed data (saved as .mat) should be loaded
%           (False; much faster).
%
%  Output:
%       A collection of cell arrays of strings, matched by row:
%           - arrayOutMicRNAName: miRBase name for the processed micro-RNA
%           - arrayOutMessRNA: HGNC symbol for the mRNA
%           - arrayOutMessRNARefSeq: RefSeq identifier for the mRNA
%       a uint32 vector, matched by row:
%           - arrayOutContPlus: the context+ score, rescaled to fit a
%               uint32 vector(*-1000)
%               
%  MATLAB Toolbox Dependencies:
%   - none
%
%  Function dependencies
%   - none
% 
% Users are encouraged to examine the TargetScan website and ensure that
%  the latest version of TargetScan (current.y 7.0) is being used 
%  appropriately:
%       http://targetscan
%
% For further information on TargetScan, please refer to:
%
% 
%
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 05/03/16
% 
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform variables/folder path pre-processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %check the extraction/pre-processing setting has been assigned
    if ~isfield(structSettings, 'processTargetScan'),
        %by default just load the data, don't re-run the extraction
        structSettings.processTargetScan = false;
    end

    %check the dB folder path has been assigned
    if ~isfield(structSettings, 'TargetScanFolder'),
        %default to the location on JCs desktop
        structSettings.TargetScanFolder = 'C:\db\targetscan_7p0\';
    end

         
    if structSettings.processTargetScan,
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Extract the data in to the output vectors
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
        
        disp(['... extracting from the TargetScan summary data file, this can ' ...
              'take some time depending on the system specifications, and it ' ...
              'requires several GB of system memory']);

        %load the TargetScan miR family file into a cell array
        filemiRContScores = fopen([structSettings.TargetScanFolder 'Summary_Counts.txt']);
        arrayInputCell = textscan(filemiRContScores,'%s','delimiter','\n');
        fclose(filemiRContScores);
        arrayCellmiRContScores = arrayInputCell{1};
        clear arrayInputCell;
        
        %examine the size of the array
        numLines = length(arrayCellmiRContScores);
        numInputRows = numLines-1;

        %only extract human (9606) interactions to save memory, so perform
        % a pre-processing run through the array and assign a Boolean flag
        arrayMicRNAIsHu = false(numInputRows,1);
        for iInt = 1:numInputRows,
            stringRow = arrayCellmiRContScores{iInt+1};
            arrayTabPos = strfind(stringRow, char(9));

            stringSpeciesID = stringRow((arrayTabPos(3)+1):(arrayTabPos(4)-1));

            if strncmp(stringSpeciesID,'9606',4) && (length(stringSpeciesID) == 4),
                arrayMicRNAIsHu(iInt) = true;
            end

        end
        %extract indices of interest from the Boolean array
        arrayMicRNAHuRows = find(arrayMicRNAIsHu)+1;
        numHuMicRNARows = length(arrayMicRNAHuRows);
        

        %create the output arrays
        arrayOutMicRNAName = cell(numHuMicRNARows,1);
        arrayOutMessRNA = cell(numHuMicRNARows,1);
        arrayOutMessRNARefSeq = cell(numHuMicRNARows,1);
        arrayOutContPlus = zeros(numHuMicRNARows,1,'uint32');

        %move through the specified output rows
        for iRow = 1:numHuMicRNARows,
            %extract the row index number
            numDataRow = arrayMicRNAHuRows(iRow);
            
            %extract the row as a string
            stringRow = arrayCellmiRContScores{numDataRow};
            %determine the delimiter (tab) positions
            arrayTabPos = strfind(stringRow, char(9));

            %mRNA name within the second column
            stringMessRNAName = stringRow((arrayTabPos(1)+1):(arrayTabPos(2)-1));
            arrayOutMessRNA{iRow} = stringMessRNAName;

            %mRNA RefSeq identifier within the first column
            stringMessRNARefSeq = stringRow(1:(arrayTabPos(1)-1));
            arrayOutMessRNARefSeq{iRow} = stringMessRNARefSeq;

            %miR name within the 14th column
            stringMicRNA = stringRow((arrayTabPos(13)+1):(arrayTabPos(14)-1));
            arrayOutMicRNAName{iRow} = stringMicRNA;
            
            %context+ score within the 15th column
            stringContPlus = stringRow((arrayTabPos(14)+1):(arrayTabPos(15)-1));
            %convert from a string to a double
            numContPlusDouble = str2double(stringContPlus);
            %multiply by 1000 and convert to a positive number to store the 
            %  decimal within an unsigned integer array
            arrayOutContPlus(iRow) = uint32(1000*-numContPlusDouble);

        end
            
        %save the output vectors for later use
        save([ structSettings.TargetScanFolder 'processedTargetScanSummaryData.mat' ], ...
                'arrayOutMicRNAName', 'arrayOutMessRNA', ...
                'arrayOutMessRNARefSeq', 'arrayOutContPlus');
        
    else
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Load the pre-processed data in to the output vectors
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
        
        disp('... loading processed TargetScan summary data for intermediate output file');
        
        %load the saved data into a temporary data structure
        structLoaded = load([ structSettings.TargetScanFolder 'processedTargetScanSummaryData.mat' ]);
        
        %extract into the required output vectors
        arrayOutMicRNAName = structLoaded.arrayOutMicRNAName;
        arrayOutMessRNA = structLoaded.arrayOutMessRNA;
        arrayOutMessRNARefSeq = structLoaded.arrayOutMessRNARefSeq;
        arrayOutContPlus = structLoaded.arrayOutContPlus;
        
    end
    
    
    
end

