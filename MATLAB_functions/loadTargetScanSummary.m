function [ arrayOutMicRNAName, arrayOutMessRNA, arrayOutMessRNARefSeq, arrayOutContPlus ] = loadTargetScanSummary( structSettings )
%LOADMIRBASE Summary of this function goes here
%   Detailed explanation goes here


%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 
% specify any fixed variables/folder paths
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 


    if ~isfield(structSettings, 'processTargetScan'),
        %by default just load the data, don't re-run the extraction
        structSettings.processTargetScan = false;
    end

    if ~isfield(structSettings, 'TargetScanFolder'),
        %default to the location on my desktop
        structSettings.TargetScanFolder = 'C:\db\targetscan_7p0\';
    end

     
%%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 
% extract data from the conserved site context+ score file
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   
    
    if structSettings.processTargetScan,
        
        disp(['... extracting from the TargetScan summary data file, this can ' ...
              'take some time depending on the system specifications, and it ' ...
              'requires several GB of system memory']);

        %load the TargetScan miR family file
        filemiRContScores = fopen([structSettings.TargetScanFolder 'Summary_Counts.txt']);
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

            stringSpeciesID = stringRow((arrayTabPos(3)+1):(arrayTabPos(4)-1));

            if strncmp(stringSpeciesID,'9606',4) && (length(stringSpeciesID) == 4),
                arrayMicRNAIsHu(iInt) = true;
            end

        end
        arrayMicRNAHuRows = find(arrayMicRNAIsHu)+1;
        numHuMicRNARows = length(arrayMicRNAHuRows);
        

        %create the output arrays
        arrayOutMicRNAName = cell(numHuMicRNARows,1);
        arrayOutMessRNA = cell(numHuMicRNARows,1);
        arrayOutMessRNARefSeq = cell(numHuMicRNARows,1);
        arrayOutContPlus = zeros(numHuMicRNARows,1,'uint32');


        for iRow = 1:numHuMicRNARows,
            numDataRow = arrayMicRNAHuRows(iRow);
            
            stringRow = arrayCellmiRContScores{numDataRow};
            arrayTabPos = strfind(stringRow, char(9));

            stringMessRNAName = stringRow((arrayTabPos(1)+1):(arrayTabPos(2)-1));
            arrayOutMessRNA{iRow} = stringMessRNAName;

            stringMessRNARefSeq = stringRow(1:(arrayTabPos(1)-1));
            arrayOutMessRNARefSeq{iRow} = stringMessRNARefSeq;

            stringMicRNA = stringRow((arrayTabPos(13)+1):(arrayTabPos(14)-1));
            arrayOutMicRNAName{iRow} = stringMicRNA;
            
            stringContPlus = stringRow((arrayTabPos(14)+1):(arrayTabPos(15)-1));
            numContPlusDouble = str2double(stringContPlus);
            %multiple by 1000 and convert to a positive number to store the 
            %  decimal within an unsigned integer array
            arrayOutContPlus(iRow) = uint32(1000*-numContPlusDouble);

        end
            
        save([ structSettings.TargetScanFolder 'processedTargetScanSummaryData.mat' ], 'arrayOutMicRNAName', 'arrayOutMessRNA', 'arrayOutMessRNARefSeq', 'arrayOutContPlus');
        
    else
        
        disp('... loading processed TargetScan summary data for intermediate output file');
        
        structLoaded = load([ structSettings.TargetScanFolder 'processedTargetScanSummaryData.mat' ]);
        
        arrayOutMicRNAName = structLoaded.arrayOutMicRNAName;
        arrayOutMessRNA = structLoaded.arrayOutMessRNA;
        arrayOutMessRNARefSeq = structLoaded.arrayOutMessRNARefSeq;
        arrayOutContPlus = structLoaded.arrayOutContPlus;
        
    end
    
    
    
end

