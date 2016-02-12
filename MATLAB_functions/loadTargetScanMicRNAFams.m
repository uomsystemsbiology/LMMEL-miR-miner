function [ arrayOutMicRNAFams, arrayOutMicRNAs ] = loadTargetScanMicRNAFams( structSettings )
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
        

        %create the output arrays
        arrayMicRNAFam = cell(numHuMicRNARows,1);
        arrayMicRNAName = cell(numHuMicRNARows,1);


        for iRow = 1:numHuMicRNARows,
            numDataRow = arrayMicRNAHuRows(iRow);

            stringRow = arrayCellmiRContScores{numDataRow};
            arrayTabPos = strfind(stringRow, char(9));

            stringMicRNA = stringRow((arrayTabPos(3)+1):(arrayTabPos(4)-1));
            arrayMicRNAName{iRow} = stringMicRNA;
            
            stringMicRNAFam = stringRow(1:(arrayTabPos(1)-1));
            arrayMicRNAFam{iRow} = stringMicRNAFam;    
        end
        arrayMicRNAFamStringLengths = zeros(length(arrayMicRNAFam),1,'uint8');
        for iMicRNAFam = 1:length(arrayMicRNAFam),
            arrayMicRNAFamStringLengths(iMicRNAFam) = length(arrayMicRNAFam{iMicRNAFam});
        end
        
        arrayOutMicRNAFams = unique(arrayMicRNAFam);
        arrayOutMicRNAs = cell(length(arrayOutMicRNAFams),1);
        for iMicRNAFam = 1:length(arrayOutMicRNAFams),
            stringMicRNAFam = arrayOutMicRNAFams{iMicRNAFam};
            arrayIndicesForFam =  strncmp(stringMicRNAFam, arrayMicRNAFam,length(stringMicRNAFam)) & ...
                                       (arrayMicRNAFamStringLengths == length(stringMicRNAFam)) ;
            arrayOutMicRNAs{iMicRNAFam} = arrayMicRNAName(arrayIndicesForFam);   
        end
                    
        save([ structSettings.TargetScanFolder 'processedTargetScanMicRNAFams.mat' ], 'arrayOutMicRNAFams', 'arrayOutMicRNAs');
        
    else
        
        disp('... loading processed TargetScan miR family data for intermediate output file');
        
        structLoaded = load([ structSettings.TargetScanFolder 'processedTargetScanMicRNAFams.mat' ]);
        
        arrayOutMicRNAFams = structLoaded.arrayOutMicRNAFams;
        arrayOutMicRNAs = structLoaded.arrayOutMicRNAs;
        
    end
    
    
    
end

