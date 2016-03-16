function [ arrayOutMicRNA, arrayOutPredMessRNAName, arrayOutPredMessRNAID, arrayOutPredIntMITGScore, arrayOutPredIntDetails ] = loadDIANAmicroT( structSettings )
%% [ arrayOutMicRNA, arrayOutPredMessRNAName, arrayOutPredMessRNAID, arrayOutPredIntMITGScore, arrayOutPredIntDetails ] = loadDIANAmicroT( structSettings )
% This function is designed to take in a structured array which specifies
%  the location of the DIANA-microT data files, and a boolean flag to 
%  control re-processing of the original data.
% NB: the miTG-score which is output from this function is rescaled to fit
%       into a uint16 format. In this case *1000 (0-1 rescaled to 0-1000).
%
%  Inputs:
%   - structSettings: a structured array which contains a number of
%           fields that control execution of different functions used
%           throughout this analysis. Fields required by this function
%           include:
%       'DIANAmicroTFolder' - a string specifying the folder path for the 
%           DIANA-microT data files
%       'processDIANAmicroT' - a Boolean flag specifying whether the
%           original data should be re-processed (True; takes a long time),
%           or whether pre-processed data (saved as .mat) should be loaded
%           (False; much faster).
%
%  Output:
%       A collection of cell arrays of strings, matched by row:
%           - arrayOutMicRNA: miRBase name for the processed micro-RNA
%           - arrayOutPredMessRNAName: HGNC symbol for the mRNA
%           - arrayOutPredIntDetails:
%       a uint32 vector, matched by row:
%           - arrayOutPredMessRNAID: ENSEMBL gene identifier for the mRNA,
%               expected prefix ENSG
%       a uint16 vector, matched by row:
%           - arrayOutPredIntMITGScore: the miTG-score rescaled for uint16
%               (i.e. 0-1 is rescaled to 0-1000)
%               
%  MATLAB Toolbox Dependencies:
%   - none
%
%  Function dependencies
%   - none
% 
% Users are encouraged to examine the DIANA Tools website and ensure that
%  the latest version of microT (currently CDS) is being used
%  appropriately:
%       http://diana.imis.athena-innovation.gr/DianaTools/
% NB: Download access for the DIANA-microT CDS data files requires the 
%       creation of an account at DIANA Tools
%
% For further information on DIANA-microT, please refer to:
%
% DIANA-microT web server v5.0: service integration into miRNA functional 
%   analysis workflows. 
%  Paraskevopoulou MD, Georgakilas G, Kostoulas N, Vlachos IS, Vergoulis T, 
%    Reczko M et al.
%  Nucleic acids research. 2013;41(Web Server issue):W169-73.
%  http://dx.doi.org/10.1093/nar/gkt393
% 
% Functional microRNA targets in protein coding sequences.
%  Reczko M, Maragkakis M, Alexiou P, Grosse I, Hatzigeorgiou AG.
%  Bioinformatics. 2012;28(6):771-6
%  http://dx.doi.org/10.1093/bioinformatics/bts043
% 
%
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 05/03/16
% 
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Perform pre-processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

    %check the extraction/pre-processing setting has been assigned
    if ~isfield(structSettings, 'processDIANAmicroT'),
        %by default just load the data, don't re-run the extraction (it
        % takes a relatively long time to extract the data into an indexed
        % format)
        structSettings.processDIANAmicroT = false;
    end
    %check the dB folder path has been assigned
    if ~isfield(structSettings, 'DIANAmicroTFolder'),
        %default to the location on JCs desktop
        structSettings.DIANAmicroTFolder = 'C:\db\dianaMicroT_CDS\';
    else
        %define the database filename accordingly for DIANA microT 4
        if strncmp(structSettings.DIANAmicroTFolder, 'C:\db\dianaMicroT_4\', length('c:\db\dianaMicroT_4\')),
            stringFileName = 'microtv4_data.csv';
        end
        %and DIANA microT CDS
        if strncmp(structSettings.DIANAmicroTFolder, 'C:\db\dianaMicroT_CDS\', length('c:\db\dianaMicroT_CDS\')),
            stringFileName = 'microT_CDS_data.csv';
        end
        
    end

 %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   % 
%% extract the DIANA-microT database
 %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %      

    if structSettings.processDIANAmicroT,

        %load the DIANA microT .csv file
        disp('Loading micro-RNA predicted targets from the DIANA-microT database');
        fileDianaMicroT = fopen([structSettings.DIANAmicroTFolder stringFileName]);
        arrayCellDIANAMicroT = textscan(fileDianaMicroT,'%s','delimiter','\n');
        fclose(fileDianaMicroT);

        %determine the total number of lines/entries
        % NB: there are a varying number of lines for each interaction, with
        % the 'details' lines beginning with "UTR" and "CDS"
        arrayDIANAMicroTUTRLines = strncmp('UTR3', arrayCellDIANAMicroT{1}, length('UTR3'));
        arrayDIANAMicroTCDSLines = strncmp('CDS', arrayCellDIANAMicroT{1}, length('CDS'));
        arrayDIANAMicroTEntryLines = ~(arrayDIANAMicroTUTRLines | arrayDIANAMicroTCDSLines);
        arrayDIANAMicroTEntryLines(1) = 0;
        arrayDianaMicroTEntryIndices = find(arrayDIANAMicroTEntryLines);

        numDIANAMicroTLines = length(arrayCellDIANAMicroT{1});
        numDIANAMicroTEntries = sum(arrayDIANAMicroTEntryLines);
        
        
        %do a quick pre-processing step to identify Hu miRs only
        disp('  limiting extraction to human micro RNAs (hsa-miR prefix)');
        arrayIsHumanMicRNAEntry = false(numDIANAMicroTEntries,1);
        %move through all DIANA-microT entries
        for iEntry = 1:numDIANAMicroTEntries,
            
            %identify the line index
            numEntryLine = arrayDianaMicroTEntryIndices(iEntry);
            if iEntry < numDIANAMicroTEntries
                numEndLine = arrayDianaMicroTEntryIndices(iEntry+1)-1;
            else
                numEndLine = numDIANAMicroTLines;
            end
            
            %extract the line
            stringInteractionLine = arrayCellDIANAMicroT{1}{numEntryLine};

            %delimit the entry with comma (,)
            arrayInteractionLine = textscan(stringInteractionLine,'%s','delimiter',',');

            %extract the miR name as the third entry
            stringMicroRNAName = arrayInteractionLine{1}{3};
            %search for bracket locations to identify the processed miR
            arrayMicroRNABracketPos = strfind(stringMicroRNAName, '(');
            if ~isempty(arrayMicroRNABracketPos),
                %extract the processed miR name
                stringProcMicroRNAName = stringMicroRNAName(1:(arrayMicroRNABracketPos(1)-1));
            end
            %check that it begins with hsa- to indicate a human miR
            if strncmp(stringProcMicroRNAName, 'hsa-', 4),
                %label the Boolean array
                arrayIsHumanMicRNAEntry(iEntry) = true;
            end
            
        end
        
        %create indices from the Boolean array with Hu entry flags
        arrayOutputEntries = find(arrayIsHumanMicRNAEntry);
        numOutputRelationships = length(arrayOutputEntries);
                       
        disp('  extracting predicted micro RNA targets');
        %create the output arrays
        arrayOutMicRNA = cell(numOutputRelationships,1);
        arrayOutPredMessRNAName = cell(numOutputRelationships,1);
        arrayOutPredMessRNAID = zeros(numOutputRelationships,1,'uint32');
        arrayOutPredIntMITGScore = zeros(numOutputRelationships,1,'uint16');
        arrayOutPredIntDetails = cell(numOutputRelationships,1);

        %create a progress indicator for tracking the database extraction
        numProgressDecile = 1;
        numProgressIndex = ceil((numProgressDecile*numOutputRelationships)/10);
        
        %move through each entry
        for iEntry = 1:numOutputRelationships,
            
            %if necessary, provide some output on the progress
            if iEntry == numProgressIndex,
                disp(['  ' num2str(numProgressDecile*10) '%']);
                numProgressDecile = numProgressDecile+1;
                numProgressIndex = ceil(numProgressDecile*numOutputRelationships/10);
            end

            %extract the specified (Hu) entry index
            numEntryLine = arrayDianaMicroTEntryIndices(arrayOutputEntries(iEntry));
            if arrayOutputEntries(iEntry) < numDIANAMicroTEntries
                numEndLine = arrayDianaMicroTEntryIndices(arrayOutputEntries(iEntry)+1)-1;
            else
                numEndLine = numDIANAMicroTLines;
            end

            %load the line
            stringInteractionLine = arrayCellDIANAMicroT{1}{numEntryLine};
            stringDetailsLine = arrayCellDIANAMicroT{1}{(numEntryLine+1):numEndLine};

            %split the line by the delimiter (,)
            arrayInteractionLine = textscan(stringInteractionLine,'%s','delimiter',',');

            %extract the mRNA name
            stringMessRNAName = arrayInteractionLine{1}{2};
            arrayMessRNANameOpenBracket = strfind(stringMessRNAName,'(');
            arrayMessRNANameCloseBracket = strfind(stringMessRNAName,')');
            stringProcMessRNAName = stringMessRNAName((arrayMessRNANameOpenBracket(1)+1):(arrayMessRNANameCloseBracket(1)-1));
            arrayOutPredMessRNAName{iEntry} = stringProcMessRNAName;
            
            %extract the mRNA identifier
            stringMessRNAID = stringMessRNAName(1:(arrayMessRNANameOpenBracket(1)-1));
            if strncmp(stringMessRNAID(1:4), 'ENSG', 4),
                %convert to double
                numMessRNAIDDouble = str2double(stringMessRNAID(4:end));
                %ensure that the identifier number will fit as uint32
                if numMessRNAIDDouble < ((2^32)-1),
                    arrayOutPredMessRNAID(iEntry) = uint32(numMessRNAIDDouble);
                else
                    if ~isnan(numMessRNAIDDouble),
                        disp('warning: ENSG ID is too big to fit as uint32');
                    end
                end
            end
            
            %extract the miR name
            stringMicroRNAName = arrayInteractionLine{1}{3};
            arrayMicroRNABracketPos = strfind(stringMicroRNAName, '(');
            stringProcMicroRNAName = stringMicroRNAName(1:(arrayMicroRNABracketPos(1)-1));
            arrayOutMicRNA{iEntry} = stringProcMicroRNAName;
            
            %output the miTG score, rescaled from 0->1, to 0->1000 and
            % saved as uint16
            arrayOutPredIntMITGScore(iEntry) = uint16(str2double(arrayInteractionLine{1}{4})*1000);
            
            %output the details for the predicted interaction
            arrayOutPredIntDetails{iEntry} = stringDetailsLine;
            
        end
            
        %save the output arrays for future use
        save( [structSettings.DIANAmicroTFolder 'processedDIANAmicroT.mat'], ...
                'arrayOutMicRNA', 'arrayOutPredMessRNAName', ...
                'arrayOutPredMessRNAID', 'arrayOutPredIntMITGScore', ...
                'arrayOutPredIntDetails' );
        
    else
        
        %load the pre-processed DIANA-microT data into an output structure
        disp('Loading processed DIANA-microT data');
        structLoaded = load([structSettings.DIANAmicroTFolder 'processedDIANAmicroT.mat']);
        
        %and then extract into the corresponding output arrays
        arrayOutMicRNA = structLoaded.arrayOutMicRNA;
        arrayOutPredMessRNAName = structLoaded.arrayOutPredMessRNAName;
        arrayOutPredMessRNAID = structLoaded.arrayOutPredMessRNAID;
        arrayOutPredIntMITGScore = structLoaded.arrayOutPredIntMITGScore;
        arrayOutPredIntDetails = structLoaded.arrayOutPredIntDetails;
        
    end
        
end

