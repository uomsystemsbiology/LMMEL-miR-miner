% This MATLAB script contains analyses for the manuscript:
%   Systems analysis of the Ludwig Melbourne melanoma cell line panel 
%    identifies miR-29b as a driver of melanoma invasiveness.
%   MC Andrews/J Cursons, DG Hurley, M Anaka, JS Cebon, A Behren, 
%    EJ Crampin
%   Genome Biology, Submitted for review Feb 2016
%   DOI: not-yet-known
%
% This script has been written to analyse data from the LM-MEL (primary
%  melanoma) cell line panel, described in:
%   The Ludwig institute for cancer research Melbourne melanoma cell line
%    panel.
%   Behren A, Anaka M, Lo PH, Vella LJ, Davis ID, Catimel J, Cardwell T, 
%    Gedye C, Hudson C, Stan R, Cebon J.
%   Pigment Cell Melanoma Res 2013, 26:597-600.
%   http://dx.doi.org/10.1111/pcmr.12097
%
% This script requires a number of functions, and utilises a number of
%  freely-accessible databases, including:
%
%   TargetScan
%
%   DIANA-microT CDS
%
%   miRTarBase
%
%   GeneOntology dB
%
%   Hoek group list of invasive/proliferative genes in melanoma
%
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%
% For further details please contact:
%
%	Dr Miles Andrews - Cancer Immunobiology Laboratory, Olivia Newton-John 
%                       Cancer Research Institute, Australia
%       miles.andrews@onjcri.org.au
%
%   Dr Joe Cursons - Systems Biology Laboratory, University of Melbourne, 
%                     Australia
%       joseph.cursons@unimelb.edu.au
%
%   Prof Edmund Crampin - Systems Biology Laboratory, University of 
%                          Melbourne, Australia
%       edmund.crampin@unimelb.edu.au
%
%   Dr Andreas Behren - Cancer Immunobiology Laboratory, Olivia Newton-John
%                        Cancer Research Institute, Australia
%       andreas.behren@onjcri.org.au
%
% This script was written by Joe Cursons - joseph.cursons@unimelb.edu.au
%
% Last Modified: 10/02/16
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Specify settings
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 

numMITGScoreThresh = 0.64;
numContPlusThresh = -0.30;

numHistBins = 100;

%specify the color (RGB) for the invasiveness clusters;
% blue [0 0 1] for low invasiveness
% red [1 0 0] for high invasiveness
arrayInvClusterColors = { [ 0 0 1 ];
                          [ 1 0 0 ] };
        
%figure layout
numFigScaleMult = 4; %scaling for on screen display
numFigOneHeight = 240; % mm
numFigOneWidth = 150; % mm
arrayFigOnePosition = [ 50 -150 numFigOneWidth*numFigScaleMult numFigOneHeight*numFigScaleMult ];
    
numFigOneAxisLabelFontSize = 10;
numFigOneAnnotationFontSize = 8;

%figure layout
numFigTwoHeight = 55; % mm
numFigTwoWidth = 160; % mm
arrayFigTwoPosition = [ 50 -150 numFigTwoWidth*numFigScaleMult numFigTwoHeight*numFigScaleMult ];

numFigXHeight = 240; % mm
numFigXWidth = 192; % mm
arrayFigXPosition = [ 50 -150 numFigXWidth*numFigScaleMult numFigXHeight*numFigScaleMult ];

numFigTwoAxisLabelFontSize = 8;
numFigTwoMarkerFontSize = 6;
  
numStatAssocSurfaceBins = 200;
numHeatMapIntensities = 2^12;

%specify the associations with miR-29b which are of particular interest
% (using the indices within the LM-MEL arrays)
                                        %   miR           mRNA
arrayMiR29bRels = { [ 478 18648 ];      %miR-29b-3p     LAMC1
                    [ 478 37081 ];      %miR-29b-3p     PPIC
                    [ 478 18690 ] };    %miR-29b-3p     LASP1


%specify the associations with other miRs which are of particular interest
% (using the indices within the LM-MEL arrays)                
                                            %   miR         mRNA
arrayOtherRelsOfInt = { [    5  18965 ];    %let-7b-5p      LIN28B
                        [  504  39576 ];    %miR-30b-5p     RUNX2
                        [  504  40213 ];    %miR-30b-5p     SERPINE1
                        [  117  17326 ];    %miR-125b-5p    IRF4
                        [  396  43286 ];    %miR-211-5p     TGFBR2
                        [  478   5571 ];    %miR-29b-3p     CDK6
                        [  478   6439 ];    %miR-29b-3p     COL4A1
                        [  478  35992 ];    %miR-29b-3p     PDGFC
                        [  283   7223 ];    %miR-17-5p      CYBRD1
                        [  509   6594 ];    %miR-30d-5p     CPE
                        [  509  17536 ];    %miR-30d-5p     JUN
                        [ 2574  12822 ];    %miR-98-5p      HIC2
                        [  244  34357 ];    %miR-146a-5p    NRAS
                        [  301  34386 ];    %miR-185-5p     NRP1
                        [  301  42200 ];    %miR-185-5p     SRPX2
                        [  354  31533 ];    %miR-199b-5p    MBP
                        [  396   1260 ];    %miR-211-5p     ANXA1
                        [  396  43085 ];    %miR-211-5p     TCF4
                        [  424   5889 ];    %miR-222-5p     CHKA
                        [  424  41829 ] };  %miR-222-5p     SOX10
                    
%sppecify the relative position for each sub-figure within Figure One                        
arrayFigOneSubPlotPositions = { [ 0.15 0.70 0.500 0.25 ];
                                [ 0.15 0.34 0.605 0.25 ];
                                [ 0.15 0.04 0.605 0.12 ];
                                [ 0.85 0.17 0.020 0.08 ] };

arrayFig2P1SubPlotPos = { [ 0.10 0.20 0.20 0.75 ];
                          [ 0.42 0.20 0.20 0.75 ];
                          [ 0.74 0.20 0.20 0.75 ] };

arrayFigXSubPlotPos = { [ 0.06 0.91 0.10 0.08 ];
                        [ 0.22 0.91 0.10 0.08 ];
                        [ 0.38 0.91 0.10 0.08 ];
                        [ 0.62 0.91 0.10 0.08 ];    %MTB: miR-125b-5p vs IRF4
                        [ 0.88 0.91 0.10 0.08 ];    %MTB
                        [ 0.14 0.67 0.10 0.08 ];    %29b CDK6
                        [ 0.46 0.67 0.10 0.08 ];    %29b COL4A1
                        [ 0.78 0.67 0.10 0.08 ];    %29b PDGFC
                        [ 0.06 0.43 0.10 0.08 ];
                        [ 0.22 0.43 0.10 0.08 ];
                        [ 0.38 0.43 0.10 0.08 ];
                        [ 0.54 0.43 0.10 0.08 ];
                        [ 0.70 0.43 0.10 0.08 ];
                        [ 0.88 0.43 0.10 0.08 ];
                        [ 0.06 0.19 0.10 0.08 ];
                        [ 0.22 0.19 0.10 0.08 ];
                        [ 0.38 0.19 0.10 0.08 ];
                        [ 0.54 0.19 0.10 0.08 ];
                        [ 0.70 0.19 0.10 0.08 ];
                        [ 0.88 0.19 0.10 0.08 ]};      

arrayInputPigmentGOTerms = { 'pigment'; 'melan' };
arrayInputEMPGOTerms = { 'epithel'; 'mesench'; 'EMT' };                       
                       
numOntologyMinGenes = 5;
numOntologyMaxGenes = 500;

%control the processing of various databases used within the script and
% specify various output folders
structSettings = struct( 'InputFolder', 'C:\wc\2015_ludwig_melanoma\data\', ...
                         'IntermediateOutputFolder', 'C:\wc\2015_ludwig_melanoma\output\', ...
                         'OutputFolder', [ stringCurrentDir ], ...
                         'flagLoadStats', true, ...
                         'flagReCalcStats', false, ...
                         'flagSaveStats', false, ...
                         'mirTarBaseFolder', 'C:\db\miRTarBase_6p1\', ...
                         'processMirTarBase', false, ...
                         'DIANAmicroTFolder', 'C:\db\dianaMicroT_CDS\', ...
                         'processDIANAmicroT', false, ...
                         'TargetScanFolder', 'C:\db\targetscan_7p0\', ...
                         'processTargetScan', false );

%initiate the script by identifying the current directory path
stringCurrentDir = cd;

%add the path to any required functions
addpath(genpath('C:\wc\code\'));

%depending upon the OS, use a different folder separater (forward vs back
% slash)
if ispc,
    strFoldSep = '\';
elseif isunix,
    strFoldSep = '/';
else
    disp(['warning: cannot determine the operating system, defaulting to ' ...
            'forward slash for the folder path separator' ] );
    strFoldSep = '/';
end

                     
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Load the Ludwig Melanoma Dataset and Perform Statistical Analysis
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 

%load in the quantitative miR/mRNA data, together with clustering 
% information 
[ structData, structClusters ] = loadLudwigData(structSettings); 
numCellLines = size(structData(1).CellLine,1);
numMicRNAs = size(structData(1).Data,1);
numMessRNAs = size(structData(2).Data,1);
%extract miR name string length arrays for subsequent searching
arrayMicRNAStringLength = zeros(length(structData(1).Target),1,'uint8');
for iMicRNA = 1:length(structData(1).Target),
    arrayMicRNAStringLength(iMicRNA) = length(structData(1).Target{iMicRNA});
end
%extract mRNA name string length arrays for subsequent searching
arrayMessRNAStringLength = zeros(length(structData(2).Target),1,'uint8');
for iMessRNA = 1:length(structData(2).Target),
    arrayMessRNAStringLength(iMessRNA) = length(structData(2).Target{iMessRNA});
end
%extract cell line string lengths for matching between miRs and mRNAs
arrayMessRNACellLineStringLengths = zeros(numCellLines,1,'int16');
for iCellLine = 1:numCellLines,
    arrayMessRNACellLineStringLengths(iCellLine) = length(structData(2).CellLine{iCellLine});
end

arrayLowInvCellLines = structClusters(2).groupMembers{1};
arrayLowInvCellLineStrLength = zeros(length(arrayLowInvCellLines),1,'uint8');
for iCellLine = 1:length(arrayLowInvCellLines),
    arrayLowInvCellLineStrLength(iCellLine) = length(arrayLowInvCellLines{iCellLine});
end

arrayHighInvCellLines = structClusters(2).groupMembers{2};
arrayHighInvCellLineStrLength = zeros(length(arrayHighInvCellLines),1,'uint8');
for iCellLine = 1:length(arrayHighInvCellLines),
    arrayHighInvCellLineStrLength(iCellLine) = length(arrayHighInvCellLines{iCellLine});
end

[ structHoekLists ] = loadHoekLists( structSettings );

%calculate measures of statistical association between pairwise
% combinations of miR/mRNA
if structSettings.flagLoadStats,
    load([structSettings.IntermediateOutputFolder 'StatAssoc.mat']);
elseif ((~structSettings.flagLoadStats) || (structSettings.flagReCalcStats)),
    if structSettings.flagSaveStats,
        [arrayMutInfo, arrayPearsCorr] = calcStatAssoc(structData);
        save([structSettings.IntermediateOutputFolder 'StatAssoc.mat'], 'arrayMutInfo', 'arrayPearsCorr');
    end
end
        
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Load various databases to enrich the analysis
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 

%etract mirTarBase 
[ arrayMTBMicRNAIDs, arrayMTBMicRNANames, arrayMTBMessRNANames, arrayMTBMessRNAEntrezIDs, arrayMTBSpecies, arrayMTBExperiment ] = loadmiRTarBase(structSettings);
%examine the micro-RNA name lengths for later indexing
arrayMTBMicRNANameLengths = zeros(length(arrayMTBMicRNANames),1,'uint16');    
for iMicRNA = 1:length(arrayMTBMicRNANames),
    arrayMTBMicRNANameLengths(iMicRNA) = length(arrayMTBMicRNANames{iMicRNA});
end    
%examine the mRNA name lengths for later indexing
arrayMTBMessRNANameLengths = zeros(length(arrayMTBMessRNANames),1,'uint16');    
for iMessRNA = 1:length(arrayMTBMessRNANames),
    arrayMTBMessRNANameLengths(iMessRNA) = length(arrayMTBMessRNANames{iMessRNA});
end    
 
%extract DIANA-microT                         
[ arrayDMTMicRNANames, arrayDMTMessRNANames, arrayDMTMessRNAIDs, arrayDMTPredIntMITGScores, arrayDMTPredIntDetails ] = loadDIANAmicroT( structSettings );
%examine the micro-RNA name lengths for later indexing
arrayDMTMicRNANameLengths = zeros(length(arrayDMTMicRNANames),1,'uint16');    
for iMicRNA = 1:length(arrayDMTMicRNANames),
    arrayDMTMicRNANameLengths(iMicRNA) = length(arrayDMTMicRNANames{iMicRNA});
end    
%examine the mRNA name lengths for later indexing
arrayDMTMessRNANameLengths = zeros(length(arrayDMTMessRNANames),1,'uint16');    
for iMessRNA = 1:length(arrayDMTMessRNANames),
    arrayDMTMessRNANameLengths(iMessRNA) = length(arrayDMTMessRNANames{iMessRNA});
end   

%extract TargetScan
[ arrayTSMicRNANames, arrayTSMessRNANames, arrayTSMessRNARefSeq, arrayTSContPlus ] = loadTargetScanSummary( structSettings );
%examine the micro-RNA name lengths for later indexing
arrayTSMicRNANameLengths = zeros(length(arrayTSMicRNANames),1,'uint16');    
for iMicRNA = 1:length(arrayTSMicRNANames),
    arrayTSMicRNANameLengths(iMicRNA) = length(arrayTSMicRNANames{iMicRNA});
end    
%examine the mRNA name lengths for later indexing
arrayTSMessRNANameLengths = zeros(length(arrayTSMessRNANames),1,'uint16');    
for iMessRNA = 1:length(arrayTSMessRNANames),
    arrayTSMessRNANameLengths(iMessRNA) = length(arrayTSMessRNANames{iMessRNA});
end    

%extract the full miR family lists from TargetScan
[ arrayTSMicRNAFams, arrayTSMicRNAsByFam ] = loadTargetScanMicRNAFams( structSettings );


%extract the full set of Gene Ontology terms for target matching
structGOSettings = struct( 'performNewExtraction', true, ...
                           'dataFolder', 'C:\db\geneontology\');
[ arrayUniqueGONums, arrayUniqueGOGenes, arrayGOMembershipMatrix ] = extractFullHumanGeneOntology( structGOSettings );

%identify the GO terms associated with the unique set of GO numbers
structGOSettings = struct( 'performNewExtraction', true, ...
                           'dataFolder', 'C:\db\geneontology\');
arrayUniqueGOTerms = matchHumanGOTermsToNums(arrayUniqueGONums, structGOSettings);


%identify GO terms which contain the specified substrings
structGOMatchingFuncSettings = struct( 'biomartFolder', 'C:\db\biomart\', ...
                                       'performNewExtraction', true );
                                   
%for pigmentation
[ ~, arrayPigGOTerm, arrayPigGONum ] = biomartEnsemblGeneToGeneOntology( arrayInputPigmentGOTerms, structGOMatchingFuncSettings );
arrayUniquePigGONums = unique(arrayPigGONum);

%and EMP
[ ~, arrayEMPGOTerm, arrayEMPGONum ] = biomartEnsemblGeneToGeneOntology( arrayInputEMPGOTerms, structGOMatchingFuncSettings );
arrayUniqueEMPGONums = unique(arrayEMPGONum);

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Perform Data Pre-Processing
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 

%the miR data contains a number of 'zero' entries (as it is from 
% Small-RNASeq) which need to be excluded
arrayNonZeroMicRNAObsIndices = find(structData(1).Data);
%before calculating summary statistics for miR ata thresholding
numNonZeroMicRNATenPercentile = prctile(structData(1).Data(arrayNonZeroMicRNAObsIndices),10);
numNonZeroMicRNANinetyPercentile = prctile(structData(1).Data(arrayNonZeroMicRNAObsIndices),90);


%calculate summary statistics for mRNA data thresholding
numMessRNATenPercentile = prctile(structData(2).Data(:),10);
numMessRNATenPercentRange = 0.1*range(structData(2).Data(:));
%use the data derived thresholds to restrict the number of miRs being
% analysed 
arrayMicRNAPassesAbundThresh = false(numMicRNAs,1);
arrayMicRNAPassesDynRangeThresh = false(numMicRNAs,1);
for iMicRNA = 1:numMicRNAs,
    if sum(structData(1).Data(iMicRNA,:) > numNonZeroMicRNATenPercentile) > 0.25*numCellLines,
        arrayMicRNAPassesAbundThresh(iMicRNA) = true;
    end
    
    if range(structData(1).Data(iMicRNA,:)) > numNonZeroMicRNANinetyPercentile,
        arrayMicRNAPassesDynRangeThresh(iMicRNA) = true;
    end
end
arrayMicRNAPassesTests = arrayMicRNAPassesDynRangeThresh & arrayMicRNAPassesAbundThresh;
    
%use the data derived thresholds to restrict the number of mRNAs being
% analysed
arrayMessRNAPassesAbundThresh = false(numMessRNAs,1);
arrayMessRNAPassesDynRangeThresh = false(numMessRNAs,1);
for iMessRNA = 1:numMessRNAs,
    if sum(structData(2).Data(iMessRNA,:) > numMessRNATenPercentile) > 0.25*numCellLines,
        arrayMessRNAPassesAbundThresh(iMessRNA) = true;
    end
    
    if range(structData(2).Data(iMessRNA,:)) > numMessRNATenPercentRange,
        arrayMessRNAPassesDynRangeThresh(iMessRNA) = true;
    end
end
arrayMessRNAPassesTests = arrayMessRNAPassesAbundThresh & arrayMessRNAPassesDynRangeThresh;

%use the 'good data' miR and mRNA list to restrict the number of pairwise
% associations being considered
arrayPairCompPassesTests = false(numMicRNAs, numMessRNAs);
for numMicRNA = find(arrayMicRNAPassesTests),
    arrayPairCompPassesTests(numMicRNA, arrayMessRNAPassesTests) = true;
end
arrayPairCompsPassTestIndices = find(arrayPairCompPassesTests);

%from the subset of pairwise associations that we are interested in, 
% calculate summary statistics for thresholding
numMutInfoThresh = prctile(arrayMutInfo(arrayPairCompsPassTestIndices), 90);
%and bin into histograms
arrayMutInfoBinEdges = linspace(min(arrayMutInfo(arrayPairCompPassesTests))-0.01,max(arrayMutInfo(arrayPairCompPassesTests))+0.01,numHistBins);
arrayMutInfoBinCentres = (arrayMutInfoBinEdges(2:end)+arrayMutInfoBinEdges(1:(end-1)))/2;
arrayMutInfoFreq = histcounts(arrayMutInfo(arrayPairCompPassesTests),arrayMutInfoBinEdges);

numPearsCorrLowThresh = prctile(arrayPearsCorr(arrayPairCompsPassTestIndices), 2.5);
numPearsCorrUpThresh = prctile(arrayPearsCorr(arrayPairCompsPassTestIndices), 97.5);


%for the pairwise associations that we are interested in (to reduce the
% search domain..), identify those with support from the various databases
arrayPairCompIsDMTPred = false(numMicRNAs, numMessRNAs);
arrayDMTAboveThreshIndices = find(arrayDMTPredIntMITGScores > uint16(numMITGScoreThresh*1000));
for iDMTRel = 1:length(arrayDMTAboveThreshIndices),
    stringDMTMicRNA = arrayDMTMicRNANames{arrayDMTAboveThreshIndices(iDMTRel)};
    stringDMTMessRNA = arrayDMTMessRNANames{arrayDMTAboveThreshIndices(iDMTRel)};
    
    arrayMicRNAMatchIndices = find( strncmp(stringDMTMicRNA, structData(1).Target, length(stringDMTMicRNA)) & ...
                                    (arrayMicRNAStringLength == length(stringDMTMicRNA) ) );
    
    arrayMessRNAMatchIndices = find( strncmp(stringDMTMessRNA, structData(2).Target, length(stringDMTMessRNA)) & ...
                                    (arrayMessRNAStringLength == length(stringDMTMessRNA) ) );
                                
    if ~isempty(arrayMicRNAMatchIndices) & ~isempty(arrayMessRNAMatchIndices),
        if (length(arrayMicRNAMatchIndices) == 1) & (length(arrayMessRNAMatchIndices) == 1),
            if arrayMicRNAPassesTests(arrayMicRNAMatchIndices) && arrayMessRNAPassesTests(arrayMessRNAMatchIndices),
                arrayPairCompIsDMTPred(arrayMicRNAMatchIndices,arrayMessRNAMatchIndices) = true;
            end
        end
        
    end
    
end
   
disp(['matching high-confidence TargetScan predicted relationships to the miR-mRNA relationship array..']);


arrayPairCompIsTSPred = false(numMicRNAs, numMessRNAs);
arrayTSAboveThreshIndices = find(arrayTSContPlus > uint16(-numContPlusThresh*1000));

numCnt = length(arrayTSAboveThreshIndices)/20;

for iTSRel = 1:length(arrayTSAboveThreshIndices),
    
    if iTSRel >= numCnt,
        numCnt = numCnt + length(arrayTSAboveThreshIndices)/20; 
        disp([num2str((numCnt/length(arrayTSAboveThreshIndices))*100) '% complete']);
    end
    stringTSMessRNA = arrayTSMessRNANames{arrayTSAboveThreshIndices(iTSRel)};
    arrayMessRNAMatchIndices = find( strncmp(stringTSMessRNA, structData(2).Target, length(stringTSMessRNA)) & ...
                                    (arrayMessRNAStringLength == length(stringTSMessRNA) ) );
    
    %match the miR back to a family
    stringTSMicRNA = arrayTSMicRNANames{arrayTSAboveThreshIndices(iTSRel)};
    arrayFamHasMatch = false(length(arrayTSMicRNAFams),1);
    for iFam = 1:length(arrayTSMicRNAFams),
        if any(strncmp(stringTSMicRNA, arrayTSMicRNAsByFam{iFam}, length(stringTSMicRNA))),
            arrayFamHasMatch(iFam) = true;
        end
    end
    numFamIndex = find(arrayFamHasMatch);
    
    if isempty(numFamIndex),
        %disp(['warning: ' stringTSMessRNA ' cannot be matched to TargetScan families']);
    else
        if length(numFamIndex) == 1,
            
            for iMemb = 1:length(arrayTSMicRNAsByFam{numFamIndex})
                stringMicRNAFamMemb = arrayTSMicRNAsByFam{numFamIndex}{iMemb};
                
                %search for the specified RNAs amongst the DIANA-microT results
                arrayMicRNAMatchIndices = find( strncmp(stringMicRNAFamMemb, structData(1).Target, length(stringMicRNAFamMemb)) & ...
                                                (arrayMicRNAStringLength == length(stringMicRNAFamMemb) ) );
                                
                if ~isempty(arrayMicRNAMatchIndices) && ~isempty(arrayMessRNAMatchIndices),
                    for iMicRNA = 1:length(arrayMicRNAMatchIndices),
                        for iMessRNA = 1:length(arrayMessRNAMatchIndices),
                            if arrayMicRNAPassesTests(arrayMicRNAMatchIndices(iMicRNA)) & arrayMessRNAPassesTests(arrayMessRNAMatchIndices(iMessRNA)),
                                arrayPairCompIsTSPred(arrayMicRNAMatchIndices(iMicRNA),arrayMessRNAMatchIndices(iMessRNA)) = true;
                            end
                        end
                    end
                else
                    %disp(['warning: cannot find high scoring TargetScan relationship: ' stringTSMicRNA ' vs. ' stringTSMessRNA ', within the LM-MEL targets']);
                end                            
                                            
            end
            
        else
            %disp(['warning: ' stringTSMessRNA ' matched to multiple TargetScan families']);
        end
    end
   
end

arrayPairCompIsMTBPred = false(numMicRNAs, numMessRNAs);
for iMTBRel = 1:length(arrayMTBMicRNANames),
    stringMTBMicRNA = arrayMTBMicRNANames{iMTBRel};
    stringMTBMessRNA = arrayMTBMessRNANames{iMTBRel};
    
    arrayMicRNAMatchIndices = find( strncmp(stringMTBMicRNA, structData(1).Target, length(stringMTBMicRNA)) & ...
                                    (arrayMicRNAStringLength == length(stringMTBMicRNA) ) );
    
    arrayMessRNAMatchIndices = find( strncmp(stringMTBMessRNA, structData(2).Target, length(stringMTBMessRNA)) & ...
                                    (arrayMessRNAStringLength == length(stringMTBMessRNA) ) );
                          
    if ~isempty(arrayMicRNAMatchIndices) & ~isempty(arrayMessRNAMatchIndices),
        if (length(arrayMicRNAMatchIndices) == 1) & (length(arrayMessRNAMatchIndices) == 1),
            if arrayMicRNAPassesTests(arrayMicRNAMatchIndices) && arrayMessRNAPassesTests(arrayMessRNAMatchIndices),
                stringExpType = arrayMTBExperiment{iMTBRel};
                        
                %check for 'strong' evidence supporting the association
                numLucRepAssPos = strfind(stringExpType,'Luciferase reporter assay');
                numQuantPCRPos = strfind(stringExpType,'qRT-PCR');
                numWestBlotPos = strfind(stringExpType,'Western blot');

                if ~isempty(numLucRepAssPos) | ~isempty(numWestBlotPos) | ~isempty(numQuantPCRPos),
                    arrayPairCompIsMTBPred(arrayMicRNAMatchIndices,arrayMessRNAMatchIndices) = true;
                end
                
            end
        end
        
    end                            
                                
end

arrayPairCompIsTSandDMTPred = arrayPairCompIsTSPred & arrayPairCompIsDMTPred;

arrayPairCompIsOnlyTSPred = (arrayPairCompIsTSPred & ~arrayPairCompIsTSandDMTPred);
arrayPairCompIsOnlyDMTPred = (arrayPairCompIsDMTPred & ~arrayPairCompIsTSandDMTPred);


%use hist3 to bin the mutual information and Pearson's correlation values
% using a regular grid
arrayStatAssocSurface = hist3(cat(2,arrayPearsCorr(arrayPairCompPassesTests), arrayMutInfo(arrayPairCompPassesTests)), [numStatAssocSurfaceBins numStatAssocSurfaceBins]);
numMaxSurfaceFreq = max(arrayStatAssocSurface(:));

arrayStatAssocSurfaceScaled = arrayStatAssocSurface/numMaxSurfaceFreq;





%identify GO terms with the specified number of component genes
arrayNumMembersPerGOTerm = sum(arrayGOMembershipMatrix,1);
arrayGOIndicesOfInterest = find( (arrayNumMembersPerGOTerm > numOntologyMinGenes) & ...
                                 (arrayNumMembersPerGOTerm < numOntologyMaxGenes) );
%extract the length of this vector for looping/iterations
numGOsOfInterest = length(arrayGOIndicesOfInterest);

%identify GO terms with the appropriate number of genes which are also
% tagged 
arrayOutputGOPigmentFlag = false(numGOsOfInterest,1);
arrayOutputGOEMPFlag = false(numGOsOfInterest,1);
for iGOTerm = 1:numGOsOfInterest,
    
    %exract the actual index from GO terms with appropriate n, and match
    % this back to the original GO number
    numGOTerm = arrayUniqueGONums(arrayGOIndicesOfInterest(iGOTerm));
    
    %compare the GO number to the 'pigmentation' list
    if any(arrayUniquePigGONums == numGOTerm),
        arrayOutputGOPigmentFlag(iGOTerm) = true;
    end
    %compare the GO number to the 'EMP' list
    if any(arrayUniqueEMPGONums == numGOTerm),
        arrayOutputGOEMPFlag(iGOTerm) = true;
    end
        
end
%extract indices from the boolean arrays
arrayOutputGOPigmentIndex = find(arrayOutputGOPigmentFlag);
arrayOutputGOEMPIndex = find(arrayOutputGOEMPFlag);


%extract gene lists for over-representation analysis
arrayOutputGOPigmentTargets = zeros(size(arrayGOMembershipMatrix,1),1,'uint8');
for iPigGOTerm = 1:length(arrayOutputGOPigmentIndex),
    arrayOutputGOPigmentTargets = arrayOutputGOPigmentTargets + uint8(arrayGOMembershipMatrix(:,arrayGOIndicesOfInterest(arrayOutputGOPigmentIndex(iPigGOTerm))));
end
arrayOutputGOPigmentTargetIndices = find(arrayOutputGOPigmentTargets > 0);

arrayLMMELMessRNAIsGOPigmentFlag = false(numMessRNAs,1);
for iMessRNA = 1:length(arrayOutputGOPigmentTargetIndices),
    stringMessRNA = arrayUniqueGOGenes{arrayOutputGOPigmentTargetIndices(iMessRNA)};
    arrayLMMELMessRNAMatch = find( strncmp(stringMessRNA,structData(2).Target,length(stringMessRNA)) & ...
                                   (arrayMessRNAStringLength == length(stringMessRNA)) );
    if ~isempty(arrayLMMELMessRNAMatch),
        arrayLMMELMessRNAIsGOPigmentFlag(arrayLMMELMessRNAMatch) = true;
    else
        
        
        if strncmp(stringMessRNA, 'PMEL', 4) && (length(stringMessRNA) == 4),
            disp([ 'hard-coded exception: ' stringMessRNA ' matched as ' ...
                   stringMessRNAToFind ' in the LM-MEL data' ]);
            stringMessRNAToFind = 'SILV';
            arrayLMMELMessRNAMatch = find( strncmp(stringMessRNAToFind,structData(2).Target,length(stringMessRNAToFind)) & ...
                                        (arrayMessRNAStringLength == length(stringMessRNAToFind)) );
            arrayLMMELMessRNAIsGOPigmentFlag(arrayLMMELMessRNAMatch) = true;
        else
            disp(['warning: ' stringMessRNA ' cannot be matched to the LM-MEL data']);
        end
        
    end
end

%extract gene lists for over-representation analysis
arrayOutputGOEMPTargets = zeros(size(arrayGOMembershipMatrix,1),1,'uint8');
for iPigGOTerm = 1:length(arrayOutputGOEMPIndex),
    arrayOutputGOEMPTargets = arrayOutputGOEMPTargets + uint8(arrayGOMembershipMatrix(:,arrayGOIndicesOfInterest(arrayOutputGOEMPIndex(iPigGOTerm))));
end
arrayOutputGOEMPTargetIndices = find(arrayOutputGOEMPTargets > 0);

arrayLMMELMessRNAIsGOEMPFlag = false(numMessRNAs,1);
for iMessRNA = 1:length(arrayOutputGOEMPTargetIndices),
    stringMessRNA = arrayUniqueGOGenes{arrayOutputGOEMPTargetIndices(iMessRNA)};
    arrayLMMELMessRNAMatch = find( strncmp(stringMessRNA,structData(2).Target,length(stringMessRNA)) & ...
                                   (arrayMessRNAStringLength == length(stringMessRNA)) );
    if ~isempty(arrayLMMELMessRNAMatch),
        arrayLMMELMessRNAIsGOEMPFlag(arrayLMMELMessRNAMatch) = true;
    else
        disp(['warning: ' stringMessRNA ' cannot be matched to the LM-MEL data']);
    end
end


arrayLMMELMessRNAIsHoekProlif = false(numMessRNAs,1);
for iMessRNA = 1:length(structHoekLists.groupMembers{1}),
    stringMessRNA = structHoekLists.groupMembers{1}{iMessRNA};
    arrayLMMELMessRNAMatch = find( strncmp(stringMessRNA,structData(2).Target,length(stringMessRNA)) & ...
                                   (arrayMessRNAStringLength == length(stringMessRNA)) );
    if ~isempty(arrayLMMELMessRNAMatch),
        arrayLMMELMessRNAIsHoekProlif(arrayLMMELMessRNAMatch) = true;
    else
        disp(['warning: ' stringMessRNA ' from the Hoek proliferative list' ...
              ' cannot be matched to the LM-MEL data']);
    end
end

arrayLMMELMessRNAIsHoekInv = false(numMessRNAs,1);
for iMessRNA = 1:length(structHoekLists.groupMembers{2}),
    stringMessRNA = structHoekLists.groupMembers{2}{iMessRNA};
    arrayLMMELMessRNAMatch = find( strncmp(stringMessRNA,structData(2).Target,length(stringMessRNA)) & ...
                                   (arrayMessRNAStringLength == length(stringMessRNA)) );
    if ~isempty(arrayLMMELMessRNAMatch),
        arrayLMMELMessRNAIsHoekInv(arrayLMMELMessRNAMatch) = true;
    else
        disp(['warning: ' stringMessRNA ' from the Hoek proliferative list' ...
              ' cannot be matched to the LM-MEL data']);
    end
end



 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Perform Data Analysis
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %

%determine the total number of miRs indexed within the TargetScan family
% list
numMicRNAsForAnalysis = 0;
for iMicRNAFam = 1:length(arrayTSMicRNAFams),
    numMicRNAsForAnalysis = numMicRNAsForAnalysis + length(arrayTSMicRNAsByFam{iMicRNAFam});
end
 
%create some output folders/arrays
arrayMicRNAFam = zeros(numMicRNAsForAnalysis, 1, 'uint16');
arrayMicRNAFamMemb = zeros(numMicRNAsForAnalysis, 1, 'uint16');

arrayLMMELIndex = zeros(numMicRNAsForAnalysis, 1, 'uint16');

arrayGoodObs = zeros(numMicRNAsForAnalysis, 1, 'uint16');
arrayPassStatTests = zeros(numMicRNAsForAnalysis, 1, 'uint16');
arrayTSTargetsAboveContPlusThresh = zeros(numMicRNAsForAnalysis, 1, 'uint16');
arrayDMTTargetsAboveMITGThresh = zeros(numMicRNAsForAnalysis, 1, 'uint16');
arrayTSTargetsAboveContPlusThreshAndPassStatTest = zeros(numMicRNAsForAnalysis, 1, 'uint16');
arrayDMTTargetsAboveMITGThreshAndPassStatTest = zeros(numMicRNAsForAnalysis, 1, 'uint16');

arrayNumGOPigTargets = zeros(numMicRNAsForAnalysis, 1, 'uint16');
arrayNumGOEMPTargets = zeros(numMicRNAsForAnalysis, 1, 'uint16');

numCounter = 1;
for iMicRNAFam = 1:length(arrayTSMicRNAFams),
    
    for iMicRNAMemb = 1:length(arrayTSMicRNAsByFam{iMicRNAFam}),
        stringMicRNA = arrayTSMicRNAsByFam{iMicRNAFam}{iMicRNAMemb};
        
        numMatchInLMMELData = find( strncmp(stringMicRNA, structData(1).Target, length(stringMicRNA)) & ...
                                   (arrayMicRNAStringLength == length(stringMicRNA)) );
        if ~isempty(numMatchInLMMELData),
            if length(numMatchInLMMELData) == 1,
                arrayLMMELIndex(numCounter) = numMatchInLMMELData;
                
                arrayMicRNAPassDataTestVector = arrayPairCompPassesTests(numMatchInLMMELData, :);
                
                arrayMicRNATSVector = arrayPairCompIsTSPred(numMatchInLMMELData, :);
                arrayMicRNADMTVector = arrayPairCompIsDMTPred(numMatchInLMMELData, :);
                arrayMicRNAMutInfoVector = arrayMutInfo(numMatchInLMMELData, :);
                arrayMicRNAPearsCorrVector = arrayPearsCorr(numMatchInLMMELData, :);
                    
                arrayMicRNAAssocPassThreshVector = ( arrayMicRNAMutInfoVector > numMutInfoThresh ) & ...
                                             ( arrayMicRNAPearsCorrVector < numPearsCorrLowThresh);
                    
                numGoodDataObs = sum(arrayMicRNAPassDataTestVector);
                numPassStatTests = sum( arrayMicRNAAssocPassThreshVector & ...
                                        arrayMicRNAPassDataTestVector );
                
                numTSTargetsAboveContPlusThresh = sum(arrayMicRNATSVector);
                numMTBTargetsAboveMITGThresh = sum(arrayMicRNADMTVector);
                
                numTSTargetsAboveContPlusThreshAndPassStatTest = sum( arrayMicRNAAssocPassThreshVector & ...
                                                                      arrayMicRNAPassDataTestVector & ...
                                                                      arrayMicRNATSVector );
                numMTBTargetsAboveMITGThreshAndPassStatTest = sum( arrayMicRNAAssocPassThreshVector & ...
                                                                   arrayMicRNAPassDataTestVector & ...
                                                                   arrayMicRNADMTVector );
                                                               
                numGOPigTargets = sum( arrayLMMELMessRNAIsGOPigmentFlag' & ...
                                       (arrayMicRNATSVector | arrayMicRNADMTVector) );
                numGOEMPTargets = sum( arrayLMMELMessRNAIsGOEMPFlag' & ...
                                       (arrayMicRNATSVector | arrayMicRNADMTVector) );
                
                arrayGoodObs(numCounter) = numGoodDataObs;
                arrayPassStatTests(numCounter) = numPassStatTests;
                arrayTSTargetsAboveContPlusThresh(numCounter) = numTSTargetsAboveContPlusThresh;
                arrayDMTTargetsAboveMITGThresh(numCounter) = numMTBTargetsAboveMITGThresh;
                arrayTSTargetsAboveContPlusThreshAndPassStatTest(numCounter) = numTSTargetsAboveContPlusThreshAndPassStatTest;
                arrayDMTTargetsAboveMITGThreshAndPassStatTest(numCounter) = numMTBTargetsAboveMITGThreshAndPassStatTest;
                arrayNumGOPigTargets(numCounter) = numGOPigTargets;
                arrayNumGOEMPTargets(numCounter) = numGOEMPTargets;
                
            else
                disp(['warning: ' stringMicRNA ' from TargetScan matched to multiple miRs from the LM-MEL data']);
            end
        end
        
        arrayMicRNAFam(numCounter) = iMicRNAFam;
        arrayMicRNAFamMemb(numCounter) = iMicRNAMemb;
        
        numCounter = numCounter+1;
    end
end

%calculate 


%% a
arrayOutputMicRNAFlag = arrayGoodObs > 0;

arrayRelEnrichmentOfActiveTSTargets = double(arrayTSTargetsAboveContPlusThreshAndPassStatTest)./double(arrayTSTargetsAboveContPlusThresh);
arrayRelEnrichmentNaNValFlag = isnan(arrayRelEnrichmentOfActiveTSTargets);
arrayRelEnrichmentOfActiveTSTargets(arrayRelEnrichmentNaNValFlag) = 0;
[~,arrayRankByTSActiveEnrichIndex] = sort(arrayRelEnrichmentOfActiveTSTargets, 'descend');

arrayRelEnrichmentOfActiveDMTTargets = double(arrayDMTTargetsAboveMITGThreshAndPassStatTest)./double(arrayDMTTargetsAboveMITGThresh);

arrayRelEnrichmentOfGOPigTerms = double(arrayNumGOPigTargets)/double(sum(arrayLMMELMessRNAIsGOPigmentFlag));
arrayRelEnrichmentOfGOEMPTerms = double(arrayNumGOEMPTargets)/double(sum(arrayLMMELMessRNAIsGOEMPFlag));

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Create the output figure - Figure 1
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %    

figOut = figure;
set(figOut, 'Position', arrayFigOnePosition);
set(figOut, 'PaperUnits', 'centimeters', 'PaperSize', [ numFigOneWidth/10, numFigOneHeight/10 ], 'PaperPosition', [ 0, 0, numFigOneWidth/10, numFigOneHeight/10 ] );


%specify the subplot axes for the Fig. 1A
subplot('Position', arrayFigOneSubPlotPositions{1});

%calculate alpha data for the surface mesh such that the low frequency data
% points are transparent
arrayAlphaData = log(arrayStatAssocSurfaceScaled+0.0001);

hold on;

%plot the mesh+contour plot
handCont = meshc(arrayStatAssocSurfaceScaled);

%control the surface transparency etc for display
handCont(1).FaceColor = 'interp';
handCont(1).EdgeColor = 'none';
handCont(1).FaceAlpha = 'interp';
handCont(1).AlphaDataMapping = 'scaled';
handCont(1).AlphaData = arrayAlphaData;
handCont(1).FaceLighting = 'phong';
handCont(1).AmbientStrength = 0.8;


%apply the x- and y- ticks and ticklabels
arraySurfacePearsCorrTicks = linspace(-1,1,5);
arraySurfaceMutInfoTicks = linspace(min(arrayMutInfo(arrayPairCompsPassTestIndices)),max(arrayMutInfo(arrayPairCompsPassTestIndices)),5);
set(gca, 'XTick', linspace(1,numStatAssocSurfaceBins,5), 'XTickLabels', sprintf('%03.2f\n',arraySurfaceMutInfoTicks));
set(gca, 'YTick', linspace(1,numStatAssocSurfaceBins,5), 'YTickLabels', sprintf('%03.2f\n',arraySurfacePearsCorrTicks));
set(gca, 'ZTick', []);
set(gca, 'FontSize', numFigOneAxisLabelFontSize);
%zlabel({'Relative';'frequency'}, 'FontSize', numFigOneAxisLabelFontSize);
xlabel('Mutual information', 'FontSize', numFigOneAxisLabelFontSize);
ylabel('Pearson''s correlation', 'FontSize', numFigOneAxisLabelFontSize);


%convert the mutual info and pearson correlation thresholds to 'binned'
% values
numMutInfoScaledVal = numStatAssocSurfaceBins*(numMutInfoThresh-min(arrayMutInfo(arrayPairCompsPassTestIndices)))/(max(arrayMutInfo(arrayPairCompsPassTestIndices))-min(arrayMutInfo(arrayPairCompsPassTestIndices)));
numPearsCorrUpThreshScaledVal = numStatAssocSurfaceBins*(numPearsCorrUpThresh-(-1))/2;
numPearsCorrLowThreshScaledVal = numStatAssocSurfaceBins*(numPearsCorrLowThresh-(-1))/2;

%draw in the threshold values for the mutual information and Pearson's
% correlation
line([0 numStatAssocSurfaceBins],[numPearsCorrUpThreshScaledVal numPearsCorrUpThreshScaledVal], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 3);
line([0 numStatAssocSurfaceBins], [numPearsCorrLowThreshScaledVal numPearsCorrLowThreshScaledVal], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 3);
line([numMutInfoScaledVal numMutInfoScaledVal], [0 numStatAssocSurfaceBins], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 3);

%create the annotation for the Pearson's correlation axis
numDoubleArrowXCent = arrayFigOneSubPlotPositions{1}(1) + arrayFigOneSubPlotPositions{1}(3) + 0.03;
numDoubleArrowYCent = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4)/2;
numDoubleArrowYLength = 0.13;
numDoubleArrowXDiff = 0.0;
numAssocDirLabel = numDoubleArrowXCent+0.05;
annotation('doublearrow', [numDoubleArrowXCent-numDoubleArrowXDiff, numDoubleArrowXCent+numDoubleArrowXDiff], [numDoubleArrowYCent-numDoubleArrowYLength, numDoubleArrowYCent+numDoubleArrowYLength]);
annotation( 'textbox', [numDoubleArrowXCent+0.02 numDoubleArrowYCent-0.05 0.1 0.1], 'String', {'Increasing';'strength';'of linear';'association'}, ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize);
annotation( 'textbox', [ numAssocDirLabel-numDoubleArrowXDiff numDoubleArrowYCent-0.17 0.1 0.1], 'String', {'Negative';'associations'}, ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize);
annotation( 'textbox', [numAssocDirLabel+numDoubleArrowXDiff numDoubleArrowYCent+0.07 0.1 0.1], 'String', {'Positive';'associations'}, ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize);
        
%create the annotation for the Mutual information axis
numArrowXStart = arrayFigOneSubPlotPositions{1}(1) + arrayFigOneSubPlotPositions{1}(3)*0.9;
numArrowYStart = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4)*1.02;
numArrowXDiff = -0.30;
numArrowYDiff = 0.0;
annotation( 'arrow', [numArrowXStart, numArrowXStart+numArrowXDiff], [numArrowYStart, numArrowYStart+numArrowYDiff] );
annotation( 'textbox', [numArrowXStart+numArrowXDiff numArrowYStart-0.03 0.30 0.10], 'String', {'Associations tend towards';'statistical independence'}, ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize);

%label the relative fraction of data points within each sextant
numFailMIPassPosPC = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) >= numPearsCorrUpThresh) & (arrayMutInfo(arrayPairCompPassesTests) < numMutInfoThresh)));
numFailMIFailPC = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) < numPearsCorrUpThresh) & (arrayPearsCorr(arrayPairCompPassesTests) > numPearsCorrLowThresh) & (arrayMutInfo(arrayPairCompPassesTests) < numMutInfoThresh)));
numFailMIPassNegPC = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) <= numPearsCorrLowThresh) & (arrayMutInfo(arrayPairCompPassesTests) < numMutInfoThresh)));
numPassMIPassPosPC = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) >= numPearsCorrUpThresh) & (arrayMutInfo(arrayPairCompPassesTests) >= numMutInfoThresh)));
numPassMIFailPC = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) < numPearsCorrUpThresh) & (arrayPearsCorr(arrayPairCompPassesTests) > numPearsCorrLowThresh) & (arrayMutInfo(arrayPairCompPassesTests) >= numMutInfoThresh)));
numPassMIPassNegPC = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) <= numPearsCorrLowThresh) & (arrayMutInfo(arrayPairCompPassesTests) >= numMutInfoThresh)));

numTotalAssoc = sum(sum((arrayPairCompPassesTests)));

%top left
numXPos = arrayFigOneSubPlotPositions{1}(1) + arrayFigOneSubPlotPositions{1}(3)*0.1;
numYPos = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4)*0.90;
annotation( 'textbox', [ numXPos numYPos 0.01 0.01], 'String', [num2str(100*numFailMIPassPosPC/numTotalAssoc, '%03.2f') '%'], ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
            'Color', 'r');
%mid left
numXPos = arrayFigOneSubPlotPositions{1}(1) + arrayFigOneSubPlotPositions{1}(3)*0.1;
numYPos = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4)*0.48;
annotation( 'textbox', [ numXPos numYPos 0.01 0.01], 'String', [num2str(100*numFailMIFailPC/numTotalAssoc, '%03.2f') '%'], ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
            'Color', 'r');
%bottom left
numXPos = arrayFigOneSubPlotPositions{1}(1) + arrayFigOneSubPlotPositions{1}(3)*0.1;
numYPos = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4)*0.05;
annotation( 'textbox', [ numXPos numYPos 0.01 0.01], 'String', [num2str(100*numFailMIPassNegPC/numTotalAssoc, '%03.2f') '%'], ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
            'Color', 'r');
%top right
numXPos = arrayFigOneSubPlotPositions{1}(1) + arrayFigOneSubPlotPositions{1}(3)*0.90;
numYPos = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4)*0.90;
annotation( 'textbox', [ numXPos numYPos 0.01 0.01], 'String', [num2str(100*numPassMIPassPosPC/numTotalAssoc, '%03.2f') '%'], ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
            'Color', 'r');
%mid right
numXPos = arrayFigOneSubPlotPositions{1}(1) + arrayFigOneSubPlotPositions{1}(3)*0.90;
numYPos = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4)*0.48;
annotation( 'textbox', [ numXPos numYPos 0.01 0.01], 'String', [num2str(100*numPassMIFailPC/numTotalAssoc, '%03.2f') '%'], ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
            'Color', 'r');
%bottom right
numXPos = arrayFigOneSubPlotPositions{1}(1) + arrayFigOneSubPlotPositions{1}(3)*0.90;
numYPos = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4)*0.05;
annotation( 'textbox', [ numXPos numYPos 0.01 0.01], 'String', [num2str(100*numPassMIPassNegPC/numTotalAssoc, '%03.2f') '%'], ...
            'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
            'Color', 'r');      
        
        
%label the sub-figure
numSubFigLabelYPos = arrayFigOneSubPlotPositions{1}(2) + arrayFigOneSubPlotPositions{1}(4);
annotation( 'textbox', [0.02 numSubFigLabelYPos 0.08 0.05], 'String', '(A)', ...
            'FontSize', numFigOneAnnotationFontSize*2, 'FontWeight', 'bold', ...
            'LineStyle', 'none' );


%specify the subplot axes for the Fig. 1B
subplot('Position',arrayFigOneSubPlotPositions{2});
arrayPlotPos = get(gca, 'Position');

hold on;

%draw in some initial points to hold the axes prior to adding the annoation
% (so that it is on the back layer)
plot(arrayMutInfo(arrayPairCompIsTSPred), arrayPearsCorr(arrayPairCompIsTSPred), 'd', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'w', 'MarkerSize', 1);
%then set the axis limits
set(gca, 'YLim', [ -1 numPearsCorrLowThresh+0.05 ]);
set(gca, 'XLim', [ numMutInfoThresh-0.1 max(arrayMutInfo(arrayPairCompPassesTests))+0.1 ]);

%insert the sub-plot label
arrayXLim = get(gca, 'XLim');
arrayYLim = get(gca, 'YLim');

arrayFig1BAnnotationXYDiff = { [ 0.10 -0.015 ];
                               [ 0.15 -0.005 ];
                               [ 0.17  0.020 ] };
                           
%draw in text arrows to highlight specific relationships
for iRel = 1:length(arrayMiR29bRels),
    numMicRNAIndex = arrayMiR29bRels{iRel}(1);
    numMessRNAIndex = arrayMiR29bRels{iRel}(2);
    
    numPearsCorr = arrayPearsCorr(numMicRNAIndex,numMessRNAIndex);
    numMutInfo = arrayMutInfo(numMicRNAIndex,numMessRNAIndex);
    
    %numXPosScaled = (numMutInfo-arrayXLim(1)/range(arrayXLim))*arrayPlotPos(3) + arrayPlotPos(1);
    numXPosScaled = ((numMutInfo-arrayXLim(1))/range(arrayXLim))*arrayPlotPos(3) + arrayPlotPos(1);
    numYPosScaled = ((numPearsCorr-arrayYLim(1))/range(arrayYLim))*arrayPlotPos(4) + arrayPlotPos(2);
    
    stringMicRNA = structData(1).Target{numMicRNAIndex};
    stringMessRNA = structData(2).Target{numMessRNAIndex};
    
    annotation('textarrow', [numXPosScaled+arrayFig1BAnnotationXYDiff{iRel}(1) numXPosScaled], [numYPosScaled+arrayFig1BAnnotationXYDiff{iRel}(2) numYPosScaled], 'String', [stringMicRNA ' vs ' stringMessRNA], ...
                'HeadLength', 1, 'HeadWidth', 1, 'FontSize', numFigOneAnnotationFontSize);
end

%re-draw the data points (over the annotation)
handleDMTDataPoints = plot(arrayMutInfo(arrayPairCompIsOnlyDMTPred), arrayPearsCorr(arrayPairCompIsOnlyDMTPred), 's', 'MarkerEdgeColor', [ 0.0 0.0 0.7 ], 'MarkerFaceColor', [ 0.3 0.3 1.0 ], 'MarkerSize', 3);
handleTSDataPoints = plot(arrayMutInfo(arrayPairCompIsOnlyTSPred), arrayPearsCorr(arrayPairCompIsOnlyTSPred), 's', 'MarkerEdgeColor', [ 0.7 0.0 0.0 ], 'MarkerFaceColor', [ 1.0 0.3 0.3 ], 'MarkerSize', 3);
handleDMTandTSDataPoints = plot(arrayMutInfo(arrayPairCompIsTSandDMTPred), arrayPearsCorr(arrayPairCompIsTSandDMTPred), 'd', 'MarkerEdgeColor', [ 0.4 0.0 0.4 ], 'MarkerFaceColor', [ 0.8 0.3 0.8 ], 'MarkerSize', 3);
handleMTBDataPoints = plot(arrayMutInfo(arrayPairCompIsMTBPred), arrayPearsCorr(arrayPairCompIsMTBPred), 'pentagram', 'MarkerEdgeColor', [ 0.0 0.7 0.0 ], 'MarkerFaceColor', 'none', 'MarkerSize', 6);

%draw in the mutual information and Pearson's correlation thresholds
plot([min(arrayMutInfo(arrayPairCompPassesTests)) max(arrayMutInfo(arrayPairCompPassesTests))],[numPearsCorrLowThresh numPearsCorrLowThresh], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 3);
plot([numMutInfoThresh numMutInfoThresh], [-1 1], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 3);

%create the figure legend
handLegend = legend( [handleMTBDataPoints handleTSDataPoints handleDMTDataPoints handleDMTandTSDataPoints ], ...
                     'miRTarBase', 'TargetScan', 'DIANA-microT', 'TargetScan & DIANA-microT', ...
                     'Location', 'southeast', 'FontSize', numFigOneAnnotationFontSize);

%move the legend to the right a bit
arrayLegendPos = get(handLegend, 'Position');
set(handLegend, 'Position', (arrayLegendPos + [0.05 0 0 0]));
% 
% %write in the relative enrichment of dB supported associations within this
% % sextant
% numBaseXPos = arrayFigOneSubPlotPositions{2}(1) + arrayFigOneSubPlotPositions{2}(3)*1.12;
% numBaseYPos = arrayFigOneSubPlotPositions{2}(2) + arrayFigOneSubPlotPositions{2}(4)*0.07;
% numYSpacer = 0.017;
% numPassMIPassNegPCWithTSandDMTSupp = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) <= numPearsCorrLowThresh) & (arrayMutInfo(arrayPairCompPassesTests) >= numMutInfoThresh) & arrayPairCompIsTSandDMTPred(arrayPairCompPassesTests)));
% annotation( 'textbox', [ numBaseXPos numBaseYPos 0.01 0.01], 'String', [num2str(100*numPassMIPassNegPCWithTSandDMTSupp/numPassMIPassNegPC, '%03.2f') '%'], ...
%             'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
%             'Color', 'r');
% numPassMIPassNegPCWithDMTSupp = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) <= numPearsCorrLowThresh) & (arrayMutInfo(arrayPairCompPassesTests) >= numMutInfoThresh) & arrayPairCompIsOnlyDMTPred(arrayPairCompPassesTests)));
% annotation( 'textbox', [ numBaseXPos numBaseYPos+numYSpacer 0.01 0.01], 'String', [num2str(100*numPassMIPassNegPCWithDMTSupp/numPassMIPassNegPC, '%03.2f') '%'], ...
%             'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
%             'Color', 'r');
% numPassMIPassNegPCWithTSSupp = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) <= numPearsCorrLowThresh) & (arrayMutInfo(arrayPairCompPassesTests) >= numMutInfoThresh) & arrayPairCompIsOnlyTSPred(arrayPairCompPassesTests)));
% annotation( 'textbox', [ numBaseXPos numBaseYPos+(2*numYSpacer) 0.01 0.01], 'String', [num2str(100*numPassMIPassNegPCWithTSSupp/numPassMIPassNegPC, '%03.2f') '%'], ...
%             'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
%             'Color', 'r');
% numPassMIPassNegPCWithMTBSupp = sum(sum((arrayPearsCorr(arrayPairCompPassesTests) <= numPearsCorrLowThresh) & (arrayMutInfo(arrayPairCompPassesTests) >= numMutInfoThresh) & arrayPairCompIsMTBPred(arrayPairCompPassesTests)));
% annotation( 'textbox', [ numBaseXPos numBaseYPos+(3*numYSpacer) 0.01 0.01], 'String', [num2str(100*numPassMIPassNegPCWithMTBSupp/numPassMIPassNegPC, '%03.2f') '%'], ...
%             'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', numFigOneAnnotationFontSize, ...
%             'Color', 'r');


hold off;

set(gca, 'FontSize', numFigOneAxisLabelFontSize);
xlabel('Mutual information', 'FontSize', numFigOneAxisLabelFontSize);
ylabel('Pearson''s correlation', 'FontSize', numFigOneAxisLabelFontSize);


%label the sub-figure
numSubFigLabelYPos = arrayFigOneSubPlotPositions{2}(2) + arrayFigOneSubPlotPositions{2}(4);
annotation( 'textbox', [0.02 numSubFigLabelYPos 0.08 0.05], 'String', '(B)', ...
            'FontSize', numFigOneAnnotationFontSize*2, 'FontWeight', 'bold', ...
            'LineStyle', 'none' );
        

%specify the subplot axes for the Fig. 1C
subplot('Position', arrayFigOneSubPlotPositions{3});

numOutputMicRNAs = 15;

% arrayOutputMicRNAIndices = arrayRankByTSActiveEnrichIndex(1:numOutputMicRNAs);

arrayCombDataForHeatMap = zeros(6,sum(arrayOutputMicRNAFlag));
arrayCombDataForHeatMap(1,:) = arrayRelEnrichmentOfActiveTSTargets(arrayRankByTSActiveEnrichIndex(1:sum(arrayOutputMicRNAFlag)));
arrayCombDataForHeatMap(2,:) = arrayRelEnrichmentOfActiveDMTTargets(arrayRankByTSActiveEnrichIndex(1:sum(arrayOutputMicRNAFlag)));

arrayCombDataForHeatMap(5,:) = arrayRelEnrichmentOfGOEMPTerms(arrayRankByTSActiveEnrichIndex(1:sum(arrayOutputMicRNAFlag)));
arrayCombDataForHeatMap(6,:) = arrayRelEnrichmentOfGOPigTerms(arrayRankByTSActiveEnrichIndex(1:sum(arrayOutputMicRNAFlag)));

structHeatMapSettings = struct('Type','data', 'Thresh', 0.55);
[ arrayCombHeatMap, arrayLUTColored, structLUTInfo ] = createHeatMapArrays( arrayCombDataForHeatMap*100, structHeatMapSettings );

arrayCombHeatMap(3:4,:,:) = 255;

image(arrayCombHeatMap(:,1:numOutputMicRNAs,:));
set(gca, 'XTick', [], 'YTick', []);
%label the individual miRs along the x-axis
for iMicRNA = 1:numOutputMicRNAs,
    numOutputIndex = arrayOutputMicRNAIndices(iMicRNA);
    stringMicRNA = arrayTSMicRNAsByFam{arrayMicRNAFam(numOutputIndex)}{arrayMicRNAFamMemb(numOutputIndex)};
    if strncmp('hsa-', stringMicRNA, length('hsa-')),
        stringMicRNAToDisp = stringMicRNA(5:end);
    else
        stringMicRNAToDisp = stringMicRNA;
    end
    text( iMicRNA, -0.1, stringMicRNAToDisp, ...
          'FontSize', 10, 'Rotation', 60, ...
          'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom' );
end
%label the rows along the y-axis
text( -1.5, 1.5, {'Predicted';'target';'enrichment'}, ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
text( numOutputMicRNAs+0.6, 1, 'TargetScan', ...
      'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
text( numOutputMicRNAs+0.6, 2, 'DIANA-microT', ...
      'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
text( -1.5, 5.5, {'Target GO';'annotation';'enrichment'}, ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
text( numOutputMicRNAs+0.6, 5, 'Pigmentation', ...
      'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
text( numOutputMicRNAs+0.6, 6, {'EMP'}, ...
      'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );

%label the sub-figure
numSubFigLabelYPos = arrayFigOneSubPlotPositions{3}(2) + arrayFigOneSubPlotPositions{3}(4)+0.06;
annotation( 'textbox', [0.02 numSubFigLabelYPos 0.08 0.05], 'String', '(C)', ...
            'FontSize', numFigOneAnnotationFontSize*2, 'FontWeight', 'bold', ...
            'LineStyle', 'none' );
        
%draw in the heat map legend
subplot('Position', arrayFigOneSubPlotPositions{4});
image(arrayLUTColored);
set(gca, 'XTick', [], 'YTick', []);
for iHMVal = 1:length(structLUTInfo),
    text(2, structLUTInfo(iHMVal).LUTPix, [structLUTInfo(iHMVal).DispValue '%'], 'FontSize', numFigOneAnnotationFontSize );
end


%save the figure as a 300 dpi PNG file
print(figOut, '-r300', '-dpng', [ structSettings.OutputFolder strFoldSep 'Fig1.png' ]);

%close the figure window
close(figOut);

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Create the output figure - Figure 2A
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %    

figOut = figure;
set(figOut, 'Position', arrayFigTwoPosition);
set(figOut, 'PaperUnits', 'centimeters', 'PaperSize', [ numFigTwoWidth/10, numFigTwoHeight/10 ], 'PaperPosition', [ 0, 0, numFigTwoWidth/10, numFigTwoHeight/10 ] );


%specify the arrow/annotation x-y spacing for the LM-MEL-9 labelling
arrayFig2P1LMMEL9AnnotationXYDiff = { [  0.01  0.03 ];
                                      [  0.05  0.07 ];
                                      [  0.03  0.03 ] };
                                
%specify the arrow/annotation x-y spacing for the LM-MEL-42 labelling
arrayFig2P1LMMEL42AnnotationXYDiff = { [  0.010  0.04 ];
                                       [  0.003  0.05 ];
                                       [  0.010  0.04 ] };
                                 
%specify the arrow/annotation x-y spacing for the LM-MEL-45 labelling
arrayFig2P1LMMEL45AnnotationXYDiff = { [ 0.09 -0.030 ];
                                       [ 0.08  0.050 ];
                                       [ 0.07 -0.010 ] };

%specify the arrow/annotation x-y spacing for the LM-MEL-77 labelling
arrayFig2P1LMMEL77AnnotationXYDiff = { [ 0.07  0.05 ];
                                       [ 0.07  0.05 ];
                                       [ 0.06  0.10 ] };
                                   
arrayXLim = [ 0 0.4 ];
arrayYLimByPlot = { [ 09.5 13.0 ];
                    [ 03.8 10.5 ];
                    [ 10.0 13.5 ] };
              
for iSpecRel = 1:length(arrayMiR29bRels),
    
    numMicRNA = arrayMiR29bRels{iSpecRel}(1);
    numMessRNA = arrayMiR29bRels{iSpecRel}(2);
    
    strMicRNA = structData(1).Target{numMicRNA};
    strMessRNA = structData(2).Target{numMessRNA};
    
    arrayMicRNAData = structData(1).Data(numMicRNA, :);
    
    %extract the mRNA data and make sure that the cell lines are
    %properly matched
    arrayMessRNAData = zeros(length(structData(1).CellLine),1,'double');
    
    %and flag the cell lines that are low/high invasive
    arrayCellLineIsLowInv = false(length(structData(1).CellLine),1);
    arrayCellLineIsHighInv = false(length(structData(1).CellLine),1);
    arrayCellLineIsLMMEL9 = false(length(structData(1).CellLine),1);
    arrayCellLineIsLMMEL42 = false(length(structData(1).CellLine),1);
    arrayCellLineIsLMMEL45 = false(length(structData(1).CellLine),1);
    arrayCellLineIsLMMEL77 = false(length(structData(1).CellLine),1);
    
    %move through each miR cell line
    for iCellLine = 1:length(structData(1).CellLine),
        stringMicRNACellLine = structData(1).CellLine{iCellLine};

        %match the cell line for the mRNA data
        numCellLineForMessRNA = find( strncmp(structData(2).CellLine,stringMicRNACellLine,length(stringMicRNACellLine)) & ...
                                      (arrayMessRNACellLineStringLengths == length(stringMicRNACellLine))  );
        if ~isempty(numCellLineForMessRNA),
            arrayMessRNAData(iCellLine) = structData(2).Data(numMessRNA,numCellLineForMessRNA);
        else
            disp(['warning: unable to find properly matching cell lines for ' stringMicRNACellLine]);
        end
        
        if any( strncmp(arrayLowInvCellLines,stringMicRNACellLine,length(stringMicRNACellLine)) & ...
                (arrayLowInvCellLineStrLength == length(stringMicRNACellLine)) ),
           arrayCellLineIsLowInv(iCellLine) = true;
        end
        
        if any( strncmp(arrayHighInvCellLines,stringMicRNACellLine,length(stringMicRNACellLine)) & ...
                (arrayHighInvCellLineStrLength == length(stringMicRNACellLine)) ),
           arrayCellLineIsHighInv(iCellLine) = true;
        end
        
        if strncmp(stringMicRNACellLine, 'LM-MEL-9', length('LM-MEL-9')) & ...
                (length(stringMicRNACellLine) == length('LM-MEL-9')),
            arrayCellLineIsLMMEL9(iCellLine) = true;
        end
        
        if strncmp(stringMicRNACellLine, 'LM-MEL-42', length('LM-MEL-42')) & ...
                (length(stringMicRNACellLine) == length('LM-MEL-42')),
            arrayCellLineIsLMMEL42(iCellLine) = true;
        end
        
        if strncmp(stringMicRNACellLine, 'LM-MEL-45', length('LM-MEL-45')) & ...
                (length(stringMicRNACellLine) == length('LM-MEL-45')),
            arrayCellLineIsLMMEL45(iCellLine) = true;
        end
        
        if strncmp(stringMicRNACellLine, 'LM-MEL-77', length('LM-MEL-77')) & ...
                (length(stringMicRNACellLine) == length('LM-MEL-77')),
            arrayCellLineIsLMMEL77(iCellLine) = true;
        end

    end
    
    numLMMEL9Index = find(arrayCellLineIsLMMEL9);
    numXForLMMEL9 = arrayMicRNAData(numLMMEL9Index);
    numYForLMMEL9 = arrayMessRNAData(numLMMEL9Index);
    
    numLMMEL42Index = find(arrayCellLineIsLMMEL42);
    numXForLMMEL42 = arrayMicRNAData(numLMMEL42Index);
    numYForLMMEL42 = arrayMessRNAData(numLMMEL42Index);
    
    numLMMEL45Index = find(arrayCellLineIsLMMEL45);
    numXForLMMEL45 = arrayMicRNAData(numLMMEL45Index);
    numYForLMMEL45 = arrayMessRNAData(numLMMEL45Index);
    
    numLMMEL77Index = find(arrayCellLineIsLMMEL77);
    numXForLMMEL77 = arrayMicRNAData(numLMMEL77Index);
    numYForLMMEL77 = arrayMessRNAData(numLMMEL77Index);
    
    subplot('Position', arrayFig2P1SubPlotPos{iSpecRel});
        
    arrayInvasiveWasNotAssayed = ~(arrayCellLineIsLowInv | arrayCellLineIsHighInv);
    
    arrayLowInvIndex = find(arrayCellLineIsLowInv);
    arrayHighInvIndex = find(arrayCellLineIsHighInv);
    arrayInvNotAssIndex = find(arrayInvasiveWasNotAssayed);
    
    hold on;
    %plot the invasiveness-unknown cell lines in grey (consistent with TCGA
    % plots)
    plot( arrayMicRNAData(arrayInvNotAssIndex), arrayMessRNAData(arrayInvNotAssIndex), 'o', ...
            'MarkerEdgeColor', [ 0.3 0.3 0.3 ], 'MarkerFaceColor', [ 0.8 0.8 0.8 ]);
    %plot the low-invasiveness cell lines in blue
    plot( arrayMicRNAData(arrayLowInvIndex), arrayMessRNAData(arrayLowInvIndex), 'o', ...
            'MarkerEdgeColor', [ 0.0 0.0 1.0 ], 'MarkerFaceColor', [ 0.6 0.6 1.0 ]);
    %plot the high-invasiveness cell lines in red
    plot( arrayMicRNAData(arrayHighInvIndex), arrayMessRNAData(arrayHighInvIndex), 'o', ...
            'MarkerEdgeColor', [ 1.0 0.0 0.0 ], 'MarkerFaceColor', [ 1.0 0.6 0.6 ]);
        
    %the LM-MEL-9 marker is behind a high invasiveness dot in the second
    % plot and this is confusing
    if iSpecRel == 2,
        %replot the LM-MEL-9 data point
        plot( arrayMicRNAData(arrayCellLineIsLMMEL9), arrayMessRNAData(arrayCellLineIsLMMEL9), 'o', ...
                'MarkerEdgeColor', [ 0.3 0.3 0.3 ], 'MarkerFaceColor', [ 0.8 0.8 0.8 ]);
    end
        
    xlabel([strMicRNA ' abundance'], 'FontSize', numFigTwoAxisLabelFontSize);
    ylabel([strMessRNA ' abundance'], 'FontSize', numFigTwoAxisLabelFontSize);
    
    set(gca, 'XLim', arrayXLim);
    set(gca, 'YLim', arrayYLimByPlot{iSpecRel});
    arrayYLim = get(gca, 'YLim');
    
    numRescaledXForLMMEL9 = arrayFig2P1SubPlotPos{iSpecRel}(1) + arrayFig2P1SubPlotPos{iSpecRel}(3)*((numXForLMMEL9-arrayXLim(1))/(arrayXLim(2)-arrayXLim(1)));
    numRescaledYForLMMEL9 = arrayFig2P1SubPlotPos{iSpecRel}(2) + arrayFig2P1SubPlotPos{iSpecRel}(4)*((numYForLMMEL9-arrayYLim(1))/(arrayYLim(2)-arrayYLim(1)));
        
    numRescaledXForLMMEL42 = arrayFig2P1SubPlotPos{iSpecRel}(1) + arrayFig2P1SubPlotPos{iSpecRel}(3)*((numXForLMMEL42-arrayXLim(1))/(arrayXLim(2)-arrayXLim(1)));
    numRescaledYForLMMEL42 = arrayFig2P1SubPlotPos{iSpecRel}(2) + arrayFig2P1SubPlotPos{iSpecRel}(4)*((numYForLMMEL42-arrayYLim(1))/(arrayYLim(2)-arrayYLim(1)));
    
    numRescaledXForLMMEL77 = arrayFig2P1SubPlotPos{iSpecRel}(1) + arrayFig2P1SubPlotPos{iSpecRel}(3)*((numXForLMMEL77-arrayXLim(1))/(arrayXLim(2)-arrayXLim(1)));
    numRescaledYForLMMEL77 = arrayFig2P1SubPlotPos{iSpecRel}(2) + arrayFig2P1SubPlotPos{iSpecRel}(4)*((numYForLMMEL77-arrayYLim(1))/(arrayYLim(2)-arrayYLim(1)));
    
    numRescaledXForLMMEL45 = arrayFig2P1SubPlotPos{iSpecRel}(1) + arrayFig2P1SubPlotPos{iSpecRel}(3)*((numXForLMMEL45-arrayXLim(1))/(arrayXLim(2)-arrayXLim(1)));
    numRescaledYForLMMEL45 = arrayFig2P1SubPlotPos{iSpecRel}(2) + arrayFig2P1SubPlotPos{iSpecRel}(4)*((numYForLMMEL45-arrayYLim(1))/(arrayYLim(2)-arrayYLim(1)));
        
    %label the LM-MEL-9 cell line
    annotation('textarrow', [numRescaledXForLMMEL9+arrayFig2P1LMMEL9AnnotationXYDiff{iSpecRel}(1) numRescaledXForLMMEL9], ...
                            [numRescaledYForLMMEL9+arrayFig2P1LMMEL9AnnotationXYDiff{iSpecRel}(2) numRescaledYForLMMEL9], ...
                            'String', ['LM-MEL-9'], 'FontSize', numFigTwoAxisLabelFontSize, ...
                            'HeadLength', 1, 'HeadWidth', 1, 'LineWidth', 1);
    %label the LM-MEL-42 cell line
    annotation('textarrow', [numRescaledXForLMMEL42+arrayFig2P1LMMEL42AnnotationXYDiff{iSpecRel}(1) numRescaledXForLMMEL42], ...
                            [numRescaledYForLMMEL42+arrayFig2P1LMMEL42AnnotationXYDiff{iSpecRel}(2) numRescaledYForLMMEL42], ...
                            'String', ['LM-MEL-42'], 'FontSize', numFigTwoAxisLabelFontSize, ...
                            'HeadLength', 1, 'HeadWidth', 1, 'LineWidth', 1);
    %label the LM-MEL-77 cell line
    annotation('textarrow', [numRescaledXForLMMEL77+arrayFig2P1LMMEL77AnnotationXYDiff{iSpecRel}(1) numRescaledXForLMMEL77], ...
                            [numRescaledYForLMMEL77+arrayFig2P1LMMEL77AnnotationXYDiff{iSpecRel}(2) numRescaledYForLMMEL77], ...
                            'String', ['LM-MEL-77'], 'FontSize', numFigTwoAxisLabelFontSize, ...
                            'HeadLength', 1, 'HeadWidth', 1, 'LineWidth', 1);
    %label the LM-MEL-45 cell line
    annotation('textarrow', [numRescaledXForLMMEL45+arrayFig2P1LMMEL45AnnotationXYDiff{iSpecRel}(1) numRescaledXForLMMEL45], ...
                            [numRescaledYForLMMEL45+arrayFig2P1LMMEL45AnnotationXYDiff{iSpecRel}(2) numRescaledYForLMMEL45], ...
                            'String', ['LM-MEL-45'], 'FontSize', numFigTwoAxisLabelFontSize, ...
                            'HeadLength', 1, 'HeadWidth', 1, 'LineWidth', 1);
        
    hold off;  
    
end

print(figOut, '-r300', '-dpng', [ structSettings.OutputFolder strFoldSep 'Fig2_Part1.png' ]);
close(figOut);

 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Create the output figure - Figure X
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %    

figOut = figure;
set(figOut, 'Position', arrayFigXPosition);
set(figOut, 'PaperUnits', 'centimeters', 'PaperSize', [ numFigXWidth/10, numFigXHeight/10 ], 'PaperPosition', [ 0, 0, numFigXWidth/10, numFigXHeight/10 ] );

arrayXLimByPlot = { [ 00.0 01.7 ];      %let-7b-5p vs LIN28B
                    [ 00.0 01.0 ];      %miR-30-5p vs RUNX2
                    [ 00.0 00.8 ];      %30b-5p vs SERPINE1
                    [ 00.0 02.0 ];      %125b-5p vs IRF4
                    [ 00.0 03.5 ];      %211-5p vs TGFBR2
                    [ 00.0 00.4 ];      %29b-3p vs CDK6
                    [ 00.0 00.4 ];      %29b-3p vs COL4A1
                    [ 00.0 00.4 ];      %29b-3p vs PDGFC
                    [ 00.0 01.0 ];      %17-5p vs CYBRD1
                    [ 00.5 05.0 ];      %30d-5p vs CPE
                    [ 00.5 05.0 ];      %30d-5p vs JUN
                    [ 00.0 01.0 ];      %98-5p vs HIC2
                    [ 00.0 05.0 ];      %146a-5p vs NRAS
                    [ 00.0 00.4 ];      %185-5p vs NRP1
                    [ 00.0 00.4 ];      %185-5p vs SRPX2
                    [ 00.0 00.7 ];      %199b-5p vs MBP
                    [ 00.0 03.3 ];      %222-5p vs ANXA1
                    [ 00.0 03.5 ];      %222-5p vs TCF4
                    [ 00.0 00.25 ];     %222-5p vs CHKA1
                    [ 00.0 00.25 ] };   %222-5p vs SOX10
                          
arrayYLimByPlot = { [ 03.9 10.5 ];      %let-7b-5p vs LIN28B
                    [ 03.9 10.5 ];      %miR-30-5p vs RUNX2
                    [ 03.9 15.5 ];      %30b-5p vs SERPINE1
                    [ 03.9 10.5 ];      %125b-5p vs IRF4
                    [ 03.9 13.0 ];      %211-5p vs TGFBR2
                    [ 04.9 13.5 ];      %29b-3p vs CDK6
                    [ 03.9 15.5 ];      %29b-3p vs COL4A1
                    [ 03.9 12.5 ];      %29b-3p vs PDGFC
                    [ 03.9 13.5 ];      %17-5p vs CYBRD1
                    [ 03.9 10.5 ];      %30d-5p vs CPE
                    [ 03.9 13.0 ];      %30d-5p vs JUN
                    [ 03.9 10.0 ];      %98-5p vs HIC2
                    [ 03.9 10.5 ];      %146a-5p vs NRAS
                    [ 03.9 12.5 ];      %185-5p vs NRP1
                    [ 03.9 10.5 ];      %185-5p vs SRPX2
                    [ 03.9 13.5 ];      %199b-5p vs MBP
                    [ 03.9 14.0 ];      %222-5p vs ANXA1
                    [ 03.9 09.0 ];      %222-5p vs TCF4
                    [ 03.9 09.0 ];      %222-5p vs CHKA1
                    [ 03.9 14.0 ] };    %222-5p vs SOX10
              
for iSpecRel = 1:length(arrayOtherRelsOfInt),
    
    numMicRNA = arrayOtherRelsOfInt{iSpecRel}(1);
    numMessRNA = arrayOtherRelsOfInt{iSpecRel}(2);
    
    strMicRNA = structData(1).Target{numMicRNA};
    strMessRNA = structData(2).Target{numMessRNA};
    
    arrayMicRNAData = structData(1).Data(numMicRNA, :);
    
    %extract the mRNA data and make sure that the cell lines are
    %properly matched
    arrayMessRNAData = zeros(length(structData(1).CellLine),1,'double');
    
    %and flag the cell lines that are low/high invasive
    arrayCellLineIsLowInv = false(length(structData(1).CellLine),1);
    arrayCellLineIsHighInv = false(length(structData(1).CellLine),1);
    
    %move through each miR cell line
    for iCellLine = 1:length(structData(1).CellLine),
        stringMicRNACellLine = structData(1).CellLine{iCellLine};

        %match the cell line for the mRNA data
        numCellLineForMessRNA = find( strncmp(structData(2).CellLine,stringMicRNACellLine,length(stringMicRNACellLine)) & ...
                                      (arrayMessRNACellLineStringLengths == length(stringMicRNACellLine))  );
        if ~isempty(numCellLineForMessRNA),
            arrayMessRNAData(iCellLine) = structData(2).Data(numMessRNA,numCellLineForMessRNA);
        else
            disp(['warning: unable to find properly matching cell lines for ' stringMicRNACellLine]);
        end
        
        if any( strncmp(arrayLowInvCellLines,stringMicRNACellLine,length(stringMicRNACellLine)) & ...
                (arrayLowInvCellLineStrLength == length(stringMicRNACellLine)) ),
           arrayCellLineIsLowInv(iCellLine) = true;
        end
        
        if any( strncmp(arrayHighInvCellLines,stringMicRNACellLine,length(stringMicRNACellLine)) & ...
                (arrayHighInvCellLineStrLength == length(stringMicRNACellLine)) ),
           arrayCellLineIsHighInv(iCellLine) = true;
        end

    end

    %plot in the specified position (note that this is not a regular MxN)
    % array
    subplot('Position', arrayFigXSubPlotPos{iSpecRel});
        
    arrayInvasiveWasNotAssayed = ~(arrayCellLineIsLowInv | arrayCellLineIsHighInv);
    
    arrayLowInvIndex = find(arrayCellLineIsLowInv);
    arrayHighInvIndex = find(arrayCellLineIsHighInv);
    arrayInvNotAssIndex = find(arrayInvasiveWasNotAssayed);
    
    hold on;
    %plot the invasiveness-unknown cell lines in grey (consistent with TCGA
    % plots)
    plot( arrayMicRNAData(arrayInvNotAssIndex), arrayMessRNAData(arrayInvNotAssIndex), 'o', ...
            'MarkerEdgeColor', [ 0.3 0.3 0.3 ], 'MarkerFaceColor', [ 0.8 0.8 0.8 ]);
    %plot the low-invasiveness cell lines in blue
    plot( arrayMicRNAData(arrayLowInvIndex), arrayMessRNAData(arrayLowInvIndex), 'o', ...
            'MarkerEdgeColor', [ 0.0 0.0 1.0 ], 'MarkerFaceColor', [ 0.6 0.6 1.0 ]);
    %plot the high-invasiveness cell lines in red
    plot( arrayMicRNAData(arrayHighInvIndex), arrayMessRNAData(arrayHighInvIndex), 'o', ...
            'MarkerEdgeColor', [ 1.0 0.0 0.0 ], 'MarkerFaceColor', [ 1.0 0.6 0.6 ]);
        
    %strip off the preceding 'hsa-' to shorten the miR display name
    if strncmp(strMicRNA, 'hsa-', length('hsa-')),
        strMicRNAToDisp = strMicRNA(5:end);
    else
        strMicRNAToDisp = strMicRNA;
    end
            
    %label the x- and y-axes
    xlabel(strMicRNAToDisp, 'FontSize', numFigTwoAxisLabelFontSize);
    ylabel(strMessRNA, 'FontSize', numFigTwoAxisLabelFontSize);
    
    set(gca, 'XLim', arrayXLimByPlot{iSpecRel});
    arrayXLim = get(gca, 'XLim');
    set(gca, 'YLim', arrayYLimByPlot{iSpecRel});
    arrayYLim = get(gca, 'YLim');

    hold off;  
    
end

print(figOut, '-r300', '-dpng', [ structSettings.OutputFolder strFoldSep 'FigX.png' ]);
close(figOut);



 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%% Create the output text files
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %    

arrayPassStats = (arrayMutInfo > numMutInfoThresh) & (arrayPearsCorr < numPearsCorrLowThresh);

arrayPassStatsAndIsInDB = arrayPassStats & ...
                          (arrayPairCompIsTSPred | arrayPairCompIsDMTPred | arrayPairCompIsMTBPred);

[arrayOutputRowIndices,arrayOutputColIndices] = find(arrayPassStatsAndIsInDB);
 

fileCSVOutput = fopen([ structSettings.OutputFolder '\AdditionalFile2.csv' ], 'w+');
fprintf(fileCSVOutput, 'miR,mRNA,probeID,r_P,MI,context+,MITG-score,miRTarBase,GO:EMP,GO:Pig,HoekInv,HoekProlif\n');
for iOutput = 1:length(arrayOutputRowIndices),
    
    numMessRNA = arrayOutputColIndices(iOutput);
    numMicRNA = arrayOutputRowIndices(iOutput);
    
    stringMessRNA = structData(2).Target{numMessRNA};
    stringMicRNA =  structData(1).Target{numMicRNA};
    
    stringProbe = structData(2).TargetID{numMessRNA};
    
    numPearsCorr = arrayPearsCorr(numMicRNA, numMessRNA);
    numMutInfo = arrayMutInfo(numMicRNA, numMessRNA);
    
    flagMessRNAIsEMPAnnot = arrayLMMELMessRNAIsGOEMPFlag(numMessRNA);
    flagMessRNAIsPigmentAnnot = arrayLMMELMessRNAIsGOPigmentFlag(numMessRNA);
    
    flagMessRNAIsHoekInv = arrayLMMELMessRNAIsHoekInv(numMessRNA);
    flagMessRNAIsHoekProlif = arrayLMMELMessRNAIsHoekProlif(numMessRNA);
    
    %output
    fprintf(fileCSVOutput, [stringMicRNA ',' stringMessRNA ',' stringProbe ',' num2str(numPearsCorr, '%04.3f') ',' num2str(numMutInfo, '%04.3f') ]);
    
    %match the miR back to a family for TargetScan
    arrayFamHasMatch = false(length(arrayTSMicRNAFams),1);
    for iFam = 1:length(arrayTSMicRNAFams),
        if any(strncmp(stringMicRNA, arrayTSMicRNAsByFam{iFam}, length(stringTSMicRNA))),
            arrayFamHasMatch(iFam) = true;
        end
    end
    numFamIndex = find(arrayFamHasMatch);
        
    arrayMessRNATSMatch = strncmp(stringMessRNA, arrayTSMessRNANames, length(stringMessRNA)) & ...
                            (arrayTSMessRNANameLengths == length(stringMessRNA));
    
    numContextPlus = 1.0;
    if isempty(numFamIndex),
        disp(['warning: ' stringMicRNA ' cannot be matched to TargetScan families']);
    else
        if length(numFamIndex) == 1,
            
            for iMemb = 1:length(arrayTSMicRNAsByFam{numFamIndex})
                stringMicRNAFamMemb = arrayTSMicRNAsByFam{numFamIndex}{iMemb};
                
                arrayTSMicRNAMatch = strncmp(stringMicRNAFamMemb, arrayTSMicRNANames, length(stringMicRNAFamMemb)) & ...
                                        (arrayTSMicRNANameLengths == length(stringMicRNAFamMemb));
                
                arrayTSMatchIndices = find(arrayMessRNATSMatch & arrayTSMicRNAMatch);
                if ~isempty(arrayTSMatchIndices),
                    if length(arrayTSMatchIndices) == 1,
                        numRelContPlus = double(arrayTSContPlus(arrayTSMatchIndices))/-1000;
                    else
                        %disp('warning: multiple TargetScan matches');
                        numRelContPlus = max(double(arrayTSContPlus(arrayTSMatchIndices))/-1000);
                    end
                    if numRelContPlus < numContextPlus,
                        numContextPlus = numRelContPlus;
                    end
                end
                                            
            end
            
        else
            disp(['warning: ' stringMicRNA ' matched to multiple TargetScan families']);
        end
    end
    
    if numContextPlus > 0,
        %it wasnt found
        fprintf(fileCSVOutput, [ ',' '-' ]);
    else
        fprintf(fileCSVOutput, [ ',' num2str(numContextPlus, '%04.3f') ]);
    end
    
    
    numMITGScore = -1;
    
    arrayMessRNADMTMatch = strncmp(stringMessRNA, arrayDMTMessRNANames, length(stringMessRNA)) & ...
                            (arrayDMTMessRNANameLengths == length(stringMessRNA));
    arrayMicRNADMTMatch = strncmp(stringMicRNA, arrayDMTMicRNANames, length(stringMicRNA)) & ...
                            (arrayDMTMicRNANameLengths == length(stringMicRNA));

    arrayDMTMatchIndices = find(arrayMessRNADMTMatch & arrayMicRNADMTMatch);
    if ~isempty(arrayDMTMatchIndices),
        if length(arrayDMTMatchIndices) == 1,
            numRelMITGScore = double(arrayDMTPredIntMITGScores(arrayDMTMatchIndices))/1000;
            
        else
            
            if strncmp('DDX60', stringMessRNA, length('DDX60')) & ...
                  strncmp('hsa-miR-340-5p', stringMicRNA, length('hsa-miR-340-5p')),
                %weird bug, need to find if this is wrong in the original
                % DIANA-microT CDS csv file? It's the first one in this
                % case anyway.. (second one is DDX60L)
                numRelMITGScore = double(arrayDMTPredIntMITGScores(arrayDMTMatchIndices(1)))/1000;
            else
                if strncmp('SCHIP1', stringMessRNA, length('SCHIP1')) & ...
                  strncmp('hsa-miR-29a-3p', stringMicRNA, length('hsa-miR-29a-3p')),
                    %weird bug, need to find if this is wrong in the original
                    % DIANA-microT CDS csv file? It's the first one in this
                    % case anyway.. (second one is IQCJ-SCHIP1, a fusion product)
                    numRelMITGScore = double(arrayDMTPredIntMITGScores(arrayDMTMatchIndices(1)))/1000;
                else
                    disp('warning: multiple DIANA-microT matches');
                    keyboard
                end
            end
        end
        if numRelMITGScore > numMITGScore,
            numMITGScore = numRelMITGScore;
        end
    end
    
    if numMITGScore < 0,
        %it wasnt found
        fprintf(fileCSVOutput, [ ',' '-' ]);
    else
        fprintf(fileCSVOutput, [ ',' num2str(numMITGScore, '%04.3f') ]);
    end
    
    arrayMessRNAMTBMatch = strncmp(stringMessRNA, arrayMTBMessRNANames, length(stringMessRNA)) & ...
                            (arrayMTBMessRNANameLengths == length(stringMessRNA));
    arrayMicRNAMTBMatch = strncmp(stringMicRNA, arrayMTBMicRNANames, length(stringMicRNA)) & ...
                            (arrayMTBMicRNANameLengths == length(stringMicRNA));
    arrayMTBMatchIndices = find(arrayMessRNAMTBMatch & arrayMicRNAMTBMatch);
                 
    stringExpEvidence = '-';
    if ~isempty(arrayMTBMatchIndices),
        if length(arrayMTBMatchIndices) == 1,
            stringExpEvidence = arrayMTBExperiment{arrayMTBMatchIndices};
        else
            stringExpEvidence = arrayMTBExperiment{arrayMTBMatchIndices(1)};
            for iEv = 2:length(arrayMTBMatchIndices),
                stringExpEvidence = [stringExpEvidence ' & ' arrayMTBExperiment{arrayMTBMatchIndices(iEv)} ];
            end
        end
    end
    fprintf(fileCSVOutput, [ ',' stringExpEvidence ]);
    
    if flagMessRNAIsEMPAnnot,
        fprintf(fileCSVOutput, ',Y' );
    else
        fprintf(fileCSVOutput, ',-' );
    end
    
    if flagMessRNAIsPigmentAnnot,
        fprintf(fileCSVOutput, ',Y' );
    else
        fprintf(fileCSVOutput, ',-' );
    end
    
    if flagMessRNAIsHoekInv,
        fprintf(fileCSVOutput, ',Y' );
    else
        fprintf(fileCSVOutput, ',-' );
    end
    
    if flagMessRNAIsHoekProlif,
        fprintf(fileCSVOutput, ',Y' );
    else
        fprintf(fileCSVOutput, ',-' );
    end
    
    %print a line feed
    fprintf(fileCSVOutput, '\n');
    
end

fclose(fileCSVOutput);
 