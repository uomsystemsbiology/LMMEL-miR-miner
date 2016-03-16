function [ arrayOutMicRNAFams, arrayOutMicRNAs ] = loadTargetScanMicRNAFams( structSettings )
%% [ arrayOutMicRNAFams, arrayOutMicRNAs ] = loadTargetScanMicRNAFams( structSettings )
% This function is designed to take in a structured array which specifies
%  the location of the TargetScan data files, and a boolean flag to 
%  control re-processing of the original data.
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
%  A collection of nested cell arrays containing strings, matched by
%       row:
%   - arrayOutMicRNAFams: a cell array containing the name of the miR
%       family (grouped by processed miR homology)
%   - arrayOutMicRNAs: a nested cell array containing all processed miRs
%       which are members of this miR family (i.e. share seed sequence
%       homology)
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
%       http://targetscan.org/
%
% For further information on TargetScan, please refer to:
% 
% Most mammalian mRNAs are conserved targets of microRNAs.
%  Friedman RC, Farh KK, Burge CB, Bartel DP.  
%  Genome research. 2009;19(1):92-105.
%  http://dx.doi.org/10.1101/gr.082701.108
% 
% MicroRNA targeting specificity in mammals: determinants beyond seed
%   pairing.
%  Grimson A, Farh KK, Johnston WK, Garrett-Engele P, Lim LP, Bartel DP.
%  Molecular cell. 2007;27(1):91-105.
%  http://dx.doi.org/10.1016/j.molcel.2007.06.017
%
% Conserved seed pairing, often flanked by adenosines, indicates that
%   thousands of human genes are microRNA targets.
%  Lewis BP, Burge CB, Bartel DP.
%  Cell. 2005;120(1):15-20.
%  http://dx.doi.org/10.1016/j.cell.2004.12.035
%
%
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 11/03/16
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
        %default to the location on my desktop
        structSettings.TargetScanFolder = 'C:\db\targetscan_7p0\';
    end

    if structSettings.processTargetScan,
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Extract the family info file into output arrays
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
        
        disp(['... extracting from the TargetScan miR family data']);

        %load the TargetScan miR family file
        filemiRContScores = fopen([structSettings.TargetScanFolder 'miR_Family_Info.txt']);
        arrayInputCell = textscan(filemiRContScores,'%s','delimiter','\n');
        fclose(filemiRContScores);
        arrayCellmiRContScores = arrayInputCell{1};
        clear arrayInputCell;
        numLines = length(arrayCellmiRContScores);
        numInputRows = numLines-1;

        %only extract human (9606) interactions to save memory
        arrayMicRNAIsHu = false(numInputRows,1);
        for iInt = 1:numInputRows,
            stringRow = arrayCellmiRContScores{iInt+1};
            arrayTabPos = strfind(stringRow, char(9));

            stringSpeciesID = stringRow((arrayTabPos(2)+1):(arrayTabPos(3)-1));

            if strncmp(stringSpeciesID,'9606',4) && (length(stringSpeciesID) == 4),
                arrayMicRNAIsHu(iInt) = true;
            end

        end
        arrayMicRNAHuRows = find(arrayMicRNAIsHu)+1;
        numHuMicRNARows = length(arrayMicRNAHuRows);
        
        %create temporary arrays for storing all miRs and their family
        arrayMicRNAFam = cell(numHuMicRNARows,1);
        arrayMicRNAName = cell(numHuMicRNARows,1);

        %move through each of the rows
        for iRow = 1:numHuMicRNARows,
            %identify the output row index
            numDataRow = arrayMicRNAHuRows(iRow);

            %extract the row as a string
            stringRow = arrayCellmiRContScores{numDataRow};
            %search for the delimiter (tab == char(9))
            arrayTabPos = strfind(stringRow, char(9));

            %extract the miR name from the fourth column
            stringMicRNA = stringRow((arrayTabPos(3)+1):(arrayTabPos(4)-1));
            arrayMicRNAName{iRow} = stringMicRNA;
            
            %extract the miR family name from the first column
            stringMicRNAFam = stringRow(1:(arrayTabPos(1)-1));
            arrayMicRNAFam{iRow} = stringMicRNAFam;    
        end
        %extract the string length of all miR family entries for later
        % indexing/string comparisons
        arrayMicRNAFamStringLengths = zeros(length(arrayMicRNAFam),1,'uint8');
        for iMicRNAFam = 1:length(arrayMicRNAFam),
            arrayMicRNAFamStringLengths(iMicRNAFam) = length(arrayMicRNAFam{iMicRNAFam});
        end
        
        %convert the list of miR families to a unique list; because the
        % input was a cell array, the output is a cell array
        arrayOutMicRNAFams = unique(arrayMicRNAFam);
        %create a nested cell array with miR members of each family
        % (matched by row of the outer cell)
        arrayOutMicRNAs = cell(length(arrayOutMicRNAFams),1);
        %move through each miR family
        for iMicRNAFam = 1:length(arrayOutMicRNAFams),
            %extract the family name as a string
            stringMicRNAFam = arrayOutMicRNAFams{iMicRNAFam};
            %identify all entries from the temporary array which match this
            % string
            arrayIndicesForFam =  strncmp(stringMicRNAFam, arrayMicRNAFam,length(stringMicRNAFam)) & ...
                                       (arrayMicRNAFamStringLengths == length(stringMicRNAFam));
            %create a nested cell array with all miRs by indexing back to
            % the temporary arrays
            arrayOutMicRNAs{iMicRNAFam} = arrayMicRNAName(arrayIndicesForFam);   
        end
                    
        %save the output vectors for later use
        save([ structSettings.TargetScanFolder 'processedTargetScanMicRNAFams.mat' ], ...
                'arrayOutMicRNAFams', 'arrayOutMicRNAs');
        
    else
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Load the pre-processed family info file into output arrays
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
        
        disp('... loading processed TargetScan miR family data for intermediate output file');
        
        %load the saved data into a temporary data structure
        structLoaded = load([ structSettings.TargetScanFolder 'processedTargetScanMicRNAFams.mat' ]);
        
        %extract into the required output vectors
        arrayOutMicRNAFams = structLoaded.arrayOutMicRNAFams;
        arrayOutMicRNAs = structLoaded.arrayOutMicRNAs;
        
    end
    
    
    
end

