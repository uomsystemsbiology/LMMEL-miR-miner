function [ arrayOutMicRNA, arrayOutPredMessRNAName, arrayOutPredMessRNAID, arrayOutPredIntMITGScore, arrayOutPredIntDetails ] = loadDIANAmicroT( structSettings )
%loadmiRBase :: A function to enrich micro-RNA lists with information from
%                   DIANA-microT
%   
% Extract information from a local copy of the miRBase database and use
%   this to enrich a cell array containing a list of micro-RNAs (as
%   strings).
%   **NB** The structSettings.DIANAmicroTFolder variable must be properly matched with the
%       local version of miRNA.dat - last tested for miRBase version 20
%   
%   inputs:
%       - arrayInMicroRNAsOfInterest: cell array of micro-RNA names to
%               enrich with miRBase
%   outputs:
%       - structMirBaseOutput: a structured array containing miRBase
%               annotation enrichment for the specified micro-RNAs of
%               interest
%
% Created by Joe Cursons, 22/10/2014
% joseph (dot) cursons (at) unimelb (dot) edu (dot) au


 %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 
%% check input settings and revert to defaults if necessary
 %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %      

    if ~isfield(structSettings, 'processDIANAmicroT'),
        %by default just load the data, don't re-run the extraction
        structSettings.processDIANAmicroT = false;
    end

    if ~isfield(structSettings, 'DIANAmicroTFolder'),
        %default to the location on my desktop
        structSettings.DIANAmicroTFolder = 'C:\db\dianaMicroT_CDS\';
    else
        if strncmp(structSettings.DIANAmicroTFolder, 'C:\db\dianaMicroT_4\', length('c:\db\dianaMicroT_4\')),
            stringFileName = 'microtv4_data.csv';
        end
        
        if strncmp(structSettings.DIANAmicroTFolder, 'C:\db\dianaMicroT_CDS\', length('c:\db\dianaMicroT_CDS\')),
            stringFileName = 'microT_CDS_data.csv';
        end
        
    end

 %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 
%% extract the DIANA-microT database
 %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %      

    if structSettings.processDIANAmicroT,

        %load the miRBase .dat file
        disp('Loading micro-RNA predicted targets from the DIANA-microT database');
        fileDianaMicroT = fopen([structSettings.DIANAmicroTFolder stringFileName]);
        arrayCellDIANAMicroT = textscan(fileDianaMicroT,'%s','delimiter','\n');
        fclose(fileDianaMicroT);

        %determine the total number of lines/entries
        % NB: there are a varying number of lines for each interaction, with
        % the 'details' lines beginning with "UTR"
        arrayDIANAMicroTUTRLines = strncmp('UTR3', arrayCellDIANAMicroT{1}, length('UTR3'));
        arrayDIANAMicroTCDSLines = strncmp('CDS', arrayCellDIANAMicroT{1}, length('CDS'));
        arrayDIANAMicroTEntryLines = ~(arrayDIANAMicroTUTRLines | arrayDIANAMicroTCDSLines);
        arrayDIANAMicroTEntryLines(1) = 0;
        arrayDianaMicroTEntryIndices = find(arrayDIANAMicroTEntryLines);

        numDIANAMicroTLines = length(arrayCellDIANAMicroT{1});
        numDIANAMicroTEntries = sum(arrayDIANAMicroTEntryLines);
        
        
        %limit to human, otherwise the output is too big
        disp('  limiting extraction to human micro RNAs (hsa-miR prefix)');
        arrayIsHumanMicRNAEntry = false(numDIANAMicroTEntries,1);
        for iEntry = 1:numDIANAMicroTEntries,
            
            numEntryLine = arrayDianaMicroTEntryIndices(iEntry);
            if iEntry < numDIANAMicroTEntries
                numEndLine = arrayDianaMicroTEntryIndices(iEntry+1)-1;
            else
                numEndLine = numDIANAMicroTLines;
            end
            
            stringInteractionLine = arrayCellDIANAMicroT{1}{numEntryLine};

            arrayInteractionLine = textscan(stringInteractionLine,'%s','delimiter',',');

            stringMicroRNAName = arrayInteractionLine{1}{3};
            arrayMicroRNABracketPos = strfind(stringMicroRNAName, '(');
            if isempty(arrayMicroRNABracketPos),
                a=1;
            else
                stringProcMicroRNAName = stringMicroRNAName(1:(arrayMicroRNABracketPos(1)-1));
            end
            
            if strncmp(stringProcMicroRNAName, 'hsa-', 4),
                arrayIsHumanMicRNAEntry(iEntry) = true;
            end
            
        end
        
        arrayOutputEntries = find(arrayIsHumanMicRNAEntry);
        numOutputRelationships = length(arrayOutputEntries);
        
        
        
        disp('  extracting predicted micro RNA targets');
        arrayOutMicRNA = cell(numOutputRelationships,1);
        arrayOutPredMessRNAName = cell(numOutputRelationships,1);
        arrayOutPredMessRNAID = zeros(numOutputRelationships,1,'uint32');
        arrayOutPredIntMITGScore = zeros(numOutputRelationships,1,'uint16');
        arrayOutPredIntDetails = cell(numOutputRelationships,1);

        %perform a quick search to determine entries that are of interest
        numProgressDecile = 1;
        numProgressIndex = ceil((numProgressDecile*numOutputRelationships)/10);
        
        for iEntry = 1:numOutputRelationships,
            
            if iEntry == numProgressIndex,
                disp(['  ' num2str(numProgressDecile*10) '%']);
                numProgressDecile = numProgressDecile+1;
                numProgressIndex = ceil(numProgressDecile*numOutputRelationships/10);
            end

            numEntryLine = arrayDianaMicroTEntryIndices(arrayOutputEntries(iEntry));
            if arrayOutputEntries(iEntry) < numDIANAMicroTEntries
                numEndLine = arrayDianaMicroTEntryIndices(arrayOutputEntries(iEntry)+1)-1;
            else
                numEndLine = numDIANAMicroTLines;
            end

            stringInteractionLine = arrayCellDIANAMicroT{1}{numEntryLine};
            stringDetailsLine = arrayCellDIANAMicroT{1}{(numEntryLine+1):numEndLine};

            arrayInteractionLine = textscan(stringInteractionLine,'%s','delimiter',',');

            stringMessRNAName = arrayInteractionLine{1}{2};
            arrayMessRNANameOpenBracket = strfind(stringMessRNAName,'(');
            arrayMessRNANameCloseBracket = strfind(stringMessRNAName,')');
            stringProcMessRNAName = stringMessRNAName((arrayMessRNANameOpenBracket(1)+1):(arrayMessRNANameCloseBracket(1)-1));
            stringMessRNAID = stringMessRNAName(1:(arrayMessRNANameOpenBracket(1)-1));
            
            
            stringMicroRNAName = arrayInteractionLine{1}{3};
            arrayMicroRNABracketPos = strfind(stringMicroRNAName, '(');
            stringProcMicroRNAName = stringMicroRNAName(1:(arrayMicroRNABracketPos(1)-1));

            arrayOutMicRNA{iEntry} = stringProcMicroRNAName;
            arrayOutPredMessRNAName{iEntry} = stringProcMessRNAName;
            if strncmp(stringMessRNAID(1:4), 'ENSG', 4),
                numMessRNAIDDouble = str2double(stringMessRNAID(4:end));
                if numMessRNAIDDouble < ((2^32)-1),
                    arrayOutPredMessRNAID(iEntry) = uint32(numMessRNAIDDouble);
                else
                    if ~isnan(numMessRNAIDDouble),
                        disp('warning: ENSG ID is too big to fit as uint32');
                    end
                end
            end
            arrayOutPredIntMITGScore(iEntry) = uint16(str2double(arrayInteractionLine{1}{4})*1000);
            arrayOutPredIntDetails{iEntry} = stringDetailsLine;
            
        end
            
        
        save([structSettings.DIANAmicroTFolder 'processedDIANAmicroT.mat'], 'arrayOutMicRNA', 'arrayOutPredMessRNAName', 'arrayOutPredMessRNAID', 'arrayOutPredIntMITGScore', 'arrayOutPredIntDetails');
        
    else
        
        disp('Loading processed DIANA-microT data');
        structLoaded = load([structSettings.DIANAmicroTFolder 'processedDIANAmicroT.mat']);
        arrayOutMicRNA = structLoaded.arrayOutMicRNA;
        arrayOutPredMessRNAName = structLoaded.arrayOutPredMessRNAName;
        arrayOutPredMessRNAID = structLoaded.arrayOutPredMessRNAID;
        arrayOutPredIntMITGScore = structLoaded.arrayOutPredIntMITGScore;
        arrayOutPredIntDetails = structLoaded.arrayOutPredIntDetails;
        
    end
    
  
    


    
end

