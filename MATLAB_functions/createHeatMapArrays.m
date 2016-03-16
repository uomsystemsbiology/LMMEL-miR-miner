function [ arrayHeatMapcoloured, arrayLUTcoloured, structLUTInfo ] = createHeatMapArrays( arrayInputData, structDataInfo )
%% [ arrayHeatMapcoloured, arrayLUTcoloured, structLUTInfo ] = createHeatMapArrays( arrayInputData, structDataInfo )
% This function is designed to take in 
% a structured array which specifies
%
%  Inputs:
%   - arrayInputData
%   - structDataInfo: a structured array which controls the execution of 
%           this function. Fields required include:
%       'Type' - a string specifying the type of data to be plotted within 
%           a heatmap; this controls the number of colour gradient
%           thresholds etc. Expected values:
%         * 'corr' - correlation; two-sided metric, specify positive and
%               negative colour-gradient switches at given thresholds
%         * 'diffExpr' - differential expression; expect log2 transformed
%               value (i.e. all positive), specify colour-gradient switch 
%               at a statistical-test metric (i.e. q-value) threshold
%         * 'mutinfo' - mutual information; one-sided metric, specify  
%               positive colour-gradient switch at given threshold
%         * 'data' - raw data values; expect expression values (i.e. no
%               negative numbers); uses the same function loop as 'mutinfo'
%               (i.e. one-sided metric, single colour switch at a
%               threshold)
%         * 'data-2thresh' - raw data values; expect expression values 
%               (i.e. no negative numbers) and one-sided metric, two colour
%               switches applied at specified thresholds)
%       'Stat' - 
%       'Thresh' - diffExpr plots: a statistical test metric (i.e. q-value) 
%                       cut-off for specifying the colour-gradient switch 
%                  mutinfo/data plots: a number value speciying the 
%                       colour-gradient switch 
%       'LowThresh' - corr/data-2thresh: a number value speciying the 
%                       colour-gradient switch for the lower threshold
%       'HighThresh' - corr/data-2thresh: a number value speciying the 
%                       colour-gradient switch for the upper threshold
%         
%  Output:
%   - structCompLists: a 1D, length 1, structured array which contains
%           fields that specify information for:
%               (1) signatures from Hoek (2006) - details below
%               (2) signatures from Widmer (2012) - details below
%           These fields include:
%               'type'
%               'groupNames'
%               'groupMembers'
%               'stringLengths'
%               
%  MATLAB Toolbox Dependencies:
%   - ?
% 
% It is recommended that users (or all scientists interested in colour
%  representation) read the following information:
%   https://www.research.ibm.com/people/l/lloydt/colour/colour.HTM
%
% This MATLAB function has been released on this repository for the 
%  manuscript which is under review at BMC Genome Biology:
%   M.C. Andrews, J. Cursons et al. (2016). Systems analysis of the Ludwig 
%       Melbourne melanoma cell line panel identifies miR-29b repression as
%        a driver of melanoma invasiveness.
%   doi: not-yet-assigned
%
% This function has also been used to generate figures for:
%   Cursons/Leuchowius et al (2015). Cell Communication and Signaling.
%       http://dx.doi.org/10.1186/s12964-015-0106-x
% 
% This function was created by Joe Cursons:
%   joseph.cursons@unimelb.edu.au
%
% Last Updated: 16/03/16
%
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%% Specify Output/Display Settings for the Heat Map colouring and LUT
%   colour-switches etc
 %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
    
    if strcmp(structDataInfo.Type,'corr'),
        numLUTs = 5;
        numStepsInLUT = [50, 10, 1, 10, 50]; %AboveHighThresh, WithinHighThresh, Zero, WithinLowThresh, BelowLowThresh

        arrayLUTRange = cell(numLUTs,2);
        arrayLUTRange{1,1} = [0 1 0]*255;       %AboveHighThresh, Maxcolour
        arrayLUTRange{1,2} = [0.4 0.6 0.4]*255; %AboveHighThresh, Mincolour
        arrayLUTRange{2,1} = [0 0.75 0.5]*255;       %WithinHighThresh, Maxcolour
        arrayLUTRange{2,2} = [0.2 0.4 0.3]*255; %WithinHighThresh, Mincolour
        arrayLUTRange{3,1} = [0 0 0]*255;       %Zero, Maxcolour
        arrayLUTRange{3,2} = [0 0 0]*255;       %Zero, Mincolour
        arrayLUTRange{4,1} = [0.4 0.3 0.2]*255; %WithinLowThresh, Maxcolour
        arrayLUTRange{4,2} = [0.75 0.5 0.0]*255;       %WithinLowThresh, Mincolour
        arrayLUTRange{5,1} = [0.6 0.4 0.4]*255; %BelowLowThresh, Maxcolour
        arrayLUTRange{5,2} = [1 0 0]*255;       %BelowLowThresh, Mincolour
        
    elseif strcmp(structDataInfo.Type,'data-2thresh'),
        
        numLUTs = 4;
        numStepsInLUT = [ 50, 50, 50, 1 ]; %AboveHighThresh, WithinHighThresh,WithinLowThresh, Zero

        arrayLUTRange = cell(numLUTs,2);
        arrayLUTRange{1,1} = [0.0 0.0 1.0]*255; %AboveThresh, Maxcolour
        arrayLUTRange{1,2} = [0.6 0.1 0.6]*255; %AboveThresh, Mincolour
        arrayLUTRange{2,1} = [1.0 0.0 0.0]*255; %WithinHighThresh, Maxcolour
        arrayLUTRange{2,2} = [0.2 0.4 0.3]*255; %WithinHighThresh, Mincolour
        arrayLUTRange{3,1} = [0.0 1.0 0.0]*255; %WithinLowThresh, Maxcolour
        arrayLUTRange{3,2} = [0.1 0.6 0.6]*255; %WithinLowThresh, Mincolour
        arrayLUTRange{4,1} = [0.0 0.0 0.0]*255; %Zero, Maxcolour
        arrayLUTRange{4,2} = [0.0 0.0 0.0]*255; %Zero, Mincolour
        
    elseif strcmp(structDataInfo.Type,'diffExpr'),
        numLUTs = 5;
        numStepsInLUT = [50, 50, 1, 50, 50]; %AboveHighThresh, WithinHighThresh, Zero, WithinLowThresh, BelowLowThresh

        arrayLUTRange = cell(numLUTs,2);
        arrayLUTRange{1,1} = [1 0 0]*255;       %AboveHighThresh, Maxcolour
        arrayLUTRange{1,2} = [0.6 0.4 0.4]*255; %AboveHighThresh, Mincolour
        arrayLUTRange{2,1} = [0.75 0.5 0]*255;       %WithinHighThresh, Maxcolour
        arrayLUTRange{2,2} = [0.4 0.3 0.2]*255; %WithinHighThresh, Mincolour
        arrayLUTRange{3,1} = [0 0 0]*255;       %Zero, Maxcolour
        arrayLUTRange{3,2} = [0 0 0]*255;       %Zero, Mincolour
        arrayLUTRange{4,1} = [0.2 0.4 0.3]*255; %WithinLowThresh, Mincolour
        arrayLUTRange{4,2} = [0 0.75 0.5]*255;       %WithinLowThresh, Maxcolour
        arrayLUTRange{5,1} = [0.4 0.6 0.4]*255; %BelowLowThresh, Mincolour
        arrayLUTRange{5,2} = [0 1 0]*255;       %BelowLowThresh, Maxcolour
        
    elseif strcmp(structDataInfo.Type,'mutinfo') || strcmp(structDataInfo.Type,'data')
        numLUTs = 3;
        numStepsInLUT = [50, 10, 1]; %AboveThresh, WithinThresh, Zero
        
        arrayLUTRange = cell(numLUTs,2);
        arrayLUTRange{1,1} = [0 1 0]*255;       %AboveThresh, Maxcolour
        arrayLUTRange{1,2} = [0.4 0.6 0.4]*255; %AboveThresh, Mincolour
        arrayLUTRange{2,1} = [0 0.75 0.5]*255;       %WithinThresh, Maxcolour
        arrayLUTRange{2,2} = [0.2 0.4 0.3]*255; %WithinThresh, Mincolour
        arrayLUTRange{3,1} = [0 0 0]*255;       %Zero, Maxcolour
        arrayLUTRange{3,2} = [0 0 0]*255;       %Zero, Mincolour
                
    end
    
    arrayNaNLUT = [0.7 0.7 0.7]*255;

% -------------------------------------------------------------------------
%% Array Management
% -------------------------------------------------------------------------

    arrayHeatMapcoloured = zeros(size(arrayInputData,1), size(arrayInputData,2), 3, 'uint8');
    arrayLUTcoloured = zeros(sum(numStepsInLUT),1, 3,'uint8');
    structLUTInfo = struct('LUTPix', {}, 'DispValue',{});
    


% -------------------------------------------------------------------------
%% Populate the Output arrays
% -------------------------------------------------------------------------
    if strcmp(structDataInfo.Type,'corr'),

        %calculate the offset around the 'zero' range
        numZeroPosOffSet = structDataInfo.HighThresh/(numStepsInLUT(2)+0.5);
        numZeroNegOffSet = structDataInfo.LowThresh/(numStepsInLUT(4)+0.5);

        %specify the correlation values for the colour LUT
        arrayCorrToColLUT = cell(numLUTs,1);
        arrayCorrToColLUT{1} = linspace(1,structDataInfo.HighThresh, numStepsInLUT(1)+1);
        arrayCorrToColLUT{2} = linspace(structDataInfo.HighThresh, numZeroPosOffSet, numStepsInLUT(2)+1);
        arrayCorrToColLUT{3} = [numZeroPosOffSet numZeroNegOffSet];
        arrayCorrToColLUT{4} = linspace(numZeroNegOffSet, structDataInfo.LowThresh, numStepsInLUT(4)+1);
        arrayCorrToColLUT{5} = linspace(structDataInfo.LowThresh, -1, numStepsInLUT(5)+1);

        %specify the colours within the LUT
        arrayColLUT = cell(numLUTs,1);
        arrayColLUT{1} = cat(2, linspace(arrayLUTRange{1,1}(1), arrayLUTRange{1,2}(1), numStepsInLUT(1))', linspace(arrayLUTRange{1,1}(2), arrayLUTRange{1,2}(2), numStepsInLUT(1))', linspace(arrayLUTRange{1,1}(3), arrayLUTRange{1,2}(3), numStepsInLUT(1))');
        arrayColLUT{2} = cat(2, linspace(arrayLUTRange{2,1}(1), arrayLUTRange{2,2}(1), numStepsInLUT(2))', linspace(arrayLUTRange{2,1}(2), arrayLUTRange{2,2}(2), numStepsInLUT(2))', linspace(arrayLUTRange{2,1}(3), arrayLUTRange{2,2}(3), numStepsInLUT(2))');
        arrayColLUT{3} = arrayLUTRange{3,1};
        arrayColLUT{4} = cat(2, linspace(arrayLUTRange{4,1}(1), arrayLUTRange{4,2}(1), numStepsInLUT(4))', linspace(arrayLUTRange{4,1}(2), arrayLUTRange{4,2}(2), numStepsInLUT(4))', linspace(arrayLUTRange{4,1}(3), arrayLUTRange{4,2}(3), numStepsInLUT(4))');
        arrayColLUT{5} = cat(2, linspace(arrayLUTRange{5,1}(1), arrayLUTRange{5,2}(1), numStepsInLUT(5))', linspace(arrayLUTRange{5,1}(2), arrayLUTRange{5,2}(2), numStepsInLUT(5))', linspace(arrayLUTRange{5,1}(3), arrayLUTRange{5,2}(3), numStepsInLUT(5))');

        %move through each LUT
        iOutputLUT = 1;
        for iLUT = 1:numLUTs,

            numSteps = numStepsInLUT(iLUT);

            %identify variable pairs with a correlation within the specified
            %range
            for iStep = 1:numSteps,

                %apply the LUT
                numMaxVal = arrayCorrToColLUT{iLUT}(iStep);
                numMinVal = arrayCorrToColLUT{iLUT}(iStep+1);
                if iLUT < 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData <= numMaxVal) & (arrayInputData > numMinVal));
                elseif iLUT == 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData <= numMaxVal) & (arrayInputData >= numMinVal));
                elseif iLUT > 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData < numMaxVal) & (arrayInputData >= numMinVal));
                end

                %and colour the output array
                for iPixel = 1:length(arrayPixelRow),
                    arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayColLUT{iLUT}(iStep,:);
                end

                %output the LUT to be displayed next to the heat map
                arrayLUTcoloured(iOutputLUT,1,:) = arrayColLUT{iLUT}(iStep,:);
                iOutputLUT = iOutputLUT+1;

            end

        end

        %specify NaN values to another colour
        arrayIsNaN = isnan(arrayInputData);
        [arrayPixelRow, arrayPixelCol] = find(arrayIsNaN);
        for iPixel = 1:length(arrayPixelRow),
            arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayNaNLUT;
        end
        
        %specify change points on the LUT for display
        structLUTInfo(1).LUTPix = 0.5;
        structLUTInfo(1).DispValue = '1';
        structLUTInfo(2).LUTPix = numStepsInLUT(1)+0.5;
        structLUTInfo(2).DispValue = num2str(structDataInfo.HighThresh, '%0.3f');
        structLUTInfo(3).LUTPix = sum(numStepsInLUT(1:3))+0.5;
        structLUTInfo(3).DispValue = '0';
        structLUTInfo(4).LUTPix = sum(numStepsInLUT(1:4))+0.5;
        structLUTInfo(4).DispValue = num2str(structDataInfo.LowThresh, '%0.3f');
        structLUTInfo(5).LUTPix = sum(numStepsInLUT)+0.5;
        structLUTInfo(5).DispValue = '-1';

    elseif strcmp(structDataInfo.Type,'diffExpr'),

        %note that the diffExpr analysis maps all of the data onto a single
        % heatmap for plotting, but then breaks this apart for display 
        % (i.e. the heatmaps for data "under thresh" and "over thresh"
        % relate to the FDR, and the diffExpr values then need to be mapped
        % back on to this
        array1DFlagFailStatThresh = find(structDataInfo.Stat(:) >= structDataInfo.Thresh);
        array1DFlagPassStatThresh = find(structDataInfo.Stat(:) < structDataInfo.Thresh);
        
        %examine the diffExpr values associated with data which
        % passed/failed the statistical threshold
        arrayDiffExprFailStats = arrayInputData(array1DFlagFailStatThresh);
        arrayFailStatsNotInf = ~isinf(arrayDiffExprFailStats);
        arrayDiffExprPassStats = arrayInputData(array1DFlagPassStatThresh);
        arrayPassStatsNotInf = ~isinf(arrayDiffExprPassStats);
        

        %specify the correlation values for the colour LUT
        arrayDiffExprToColLUT = cell(numLUTs,1);
        if ~isempty(arrayPassStatsNotInf),
            arrayDiffExprToColLUT{1} = linspace(max(arrayDiffExprPassStats(arrayPassStatsNotInf)),1E-6, numStepsInLUT(1)+1);
        else
            arrayDiffExprToColLUT{1} = linspace(1,1E-6, numStepsInLUT(1)+1);
        end
        if ~isempty(arrayDiffExprFailStats),
            arrayDiffExprToColLUT{2} = linspace(max(arrayDiffExprFailStats(arrayFailStatsNotInf)),1E-6, numStepsInLUT(2)+1);
        else
            arrayDiffExprToColLUT{2} = linspace(2E-6,1E-6, numStepsInLUT(2)+1);
        end
        arrayDiffExprToColLUT{3} = [1E-6 -1E-6];
        if ~isempty(arrayDiffExprFailStats),
            arrayDiffExprToColLUT{4} = linspace(-1E-6, min(arrayDiffExprFailStats(arrayFailStatsNotInf)), numStepsInLUT(4)+1);
        else
            arrayDiffExprToColLUT{4} = linspace(-1E-6, -2E-6, numStepsInLUT(4)+1);
        end
        if ~isempty(arrayPassStatsNotInf),
            arrayDiffExprToColLUT{5} = linspace(-1E-6, min(arrayDiffExprPassStats(arrayPassStatsNotInf)), numStepsInLUT(5)+1);
        else
            arrayDiffExprToColLUT{5} = linspace(-1E-6, -1, numStepsInLUT(5)+1);
        end

        %specify the colours within the LUT
        arrayColLUT = cell(numLUTs,1);
        arrayColLUT{1} = cat(2, linspace(arrayLUTRange{1,1}(1), arrayLUTRange{1,2}(1), numStepsInLUT(1))', linspace(arrayLUTRange{1,1}(2), arrayLUTRange{1,2}(2), numStepsInLUT(1))', linspace(arrayLUTRange{1,1}(3), arrayLUTRange{1,2}(3), numStepsInLUT(1))');
        arrayColLUT{2} = cat(2, linspace(arrayLUTRange{2,1}(1), arrayLUTRange{2,2}(1), numStepsInLUT(2))', linspace(arrayLUTRange{2,1}(2), arrayLUTRange{2,2}(2), numStepsInLUT(2))', linspace(arrayLUTRange{2,1}(3), arrayLUTRange{2,2}(3), numStepsInLUT(2))');
        arrayColLUT{3} = arrayLUTRange{3,1};
        arrayColLUT{4} = cat(2, linspace(arrayLUTRange{4,1}(1), arrayLUTRange{4,2}(1), numStepsInLUT(4))', linspace(arrayLUTRange{4,1}(2), arrayLUTRange{4,2}(2), numStepsInLUT(4))', linspace(arrayLUTRange{4,1}(3), arrayLUTRange{4,2}(3), numStepsInLUT(4))');
        arrayColLUT{5} = cat(2, linspace(arrayLUTRange{5,1}(1), arrayLUTRange{5,2}(1), numStepsInLUT(5))', linspace(arrayLUTRange{5,1}(2), arrayLUTRange{5,2}(2), numStepsInLUT(5))', linspace(arrayLUTRange{5,1}(3), arrayLUTRange{5,2}(3), numStepsInLUT(5))');

        %plot all of the data points which fail the statistical threshold
        arrayFlagFailStatThresh = (structDataInfo.Stat >= structDataInfo.Thresh);
        arrayFlagPassStatThresh = (structDataInfo.Stat < structDataInfo.Thresh);

        
        %move through each LUT
        iOutputLUT = 1;
        for iLUT = 1:numLUTs,

            numSteps = numStepsInLUT(iLUT);

            %identify output positions with diffExpr within the specified
            % range
            for iStep = 1:numSteps,

                %apply the LUT
                numMaxVal = arrayDiffExprToColLUT{iLUT}(iStep);
                numMinVal = arrayDiffExprToColLUT{iLUT}(iStep+1);
                if iLUT < 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData <= numMaxVal) & (arrayInputData > numMinVal));
                elseif iLUT == 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData <= numMaxVal) & (arrayInputData >= numMinVal));
                elseif iLUT > 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData < numMaxVal) & (arrayInputData >= numMinVal));
                end

                %and colour the output array
                for iPixel = 1:length(arrayPixelRow),
                    if (arrayFlagFailStatThresh(arrayPixelRow(iPixel), arrayPixelCol(iPixel)) && ((iLUT == 2) || (iLUT == 4))),
                        arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayColLUT{iLUT}(iStep,:);
                    elseif (arrayFlagPassStatThresh(arrayPixelRow(iPixel), arrayPixelCol(iPixel)) && ((iLUT == 1) || (iLUT == 5))),
                        arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayColLUT{iLUT}(iStep,:);
                    else
                        %do nothing
                    end
                end

                %output the LUT to be displayed next to the heat map
                arrayLUTcoloured(iOutputLUT,1,:) = arrayColLUT{iLUT}(iStep,:);
                iOutputLUT = iOutputLUT+1;

            end

        end

        %specify NaN values to another colour
        arrayIsNaN = isnan(arrayInputData);
        [arrayPixelRow, arrayPixelCol] = find(arrayIsNaN);
        for iPixel = 1:length(arrayPixelRow),
            arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayNaNLUT;
        end
        
        %specify change points on the LUT for display
        structLUTInfo(1).LUTPix = 0.5;
        structLUTInfo(1).DispValue = [ num2str(max(arrayDiffExprPassStats(arrayPassStatsNotInf)), '%0.3f') '; < {\tau}_{FDR}'];
        structLUTInfo(2).LUTPix = numStepsInLUT(1)+0.5;
        structLUTInfo(2).DispValue = [ num2str(max(arrayDiffExprFailStats(arrayFailStatsNotInf)), '%0.3f') '; >= {\tau}_{FDR}'];
        structLUTInfo(3).LUTPix = sum(numStepsInLUT(1:3))+0.5;
        structLUTInfo(3).DispValue = '0';
        structLUTInfo(4).LUTPix = sum(numStepsInLUT(1:4))+0.5;
        structLUTInfo(4).DispValue = [ num2str(max(arrayDiffExprFailStats(arrayFailStatsNotInf)), '%0.3f') '; >= {\tau}_{FDR}'];
        structLUTInfo(5).LUTPix = sum(numStepsInLUT)+0.5;
        structLUTInfo(5).DispValue = [ num2str(min(arrayDiffExprPassStats(arrayPassStatsNotInf)), '%0.3f') '; < {\tau}_{FDR}'];

    elseif strcmp(structDataInfo.Type,'data-2thresh'),

        numZeroOffSet = max(arrayInputData(:))/10000;
        
        %specify the values for the colour LUT
        arrayDataToColLUT = cell(numLUTs,1);
        arrayDataToColLUT{1} = linspace(max(arrayInputData(:)),structDataInfo.HighThresh, numStepsInLUT(1)+1);
        arrayDataToColLUT{2} = linspace(structDataInfo.HighThresh, structDataInfo.LowThresh, numStepsInLUT(2)+1);
        arrayDataToColLUT{3} = linspace(structDataInfo.LowThresh, numZeroOffSet, numStepsInLUT(3)+1);
        arrayDataToColLUT{4} = [numZeroOffSet 0];
        
        %specify the colours within the LUT
        arrayColLUT = cell(numLUTs,1);
        arrayColLUT{1} = cat(2, linspace(arrayLUTRange{1,1}(1), arrayLUTRange{1,2}(1), numStepsInLUT(1))', linspace(arrayLUTRange{1,1}(2), arrayLUTRange{1,2}(2), numStepsInLUT(1))', linspace(arrayLUTRange{1,1}(3), arrayLUTRange{1,2}(3), numStepsInLUT(1))');
        arrayColLUT{2} = cat(2, linspace(arrayLUTRange{2,1}(1), arrayLUTRange{2,2}(1), numStepsInLUT(2))', linspace(arrayLUTRange{2,1}(2), arrayLUTRange{2,2}(2), numStepsInLUT(2))', linspace(arrayLUTRange{2,1}(3), arrayLUTRange{2,2}(3), numStepsInLUT(2))');
        arrayColLUT{3} = cat(2, linspace(arrayLUTRange{3,1}(1), arrayLUTRange{3,2}(1), numStepsInLUT(3))', linspace(arrayLUTRange{3,1}(2), arrayLUTRange{3,2}(2), numStepsInLUT(3))', linspace(arrayLUTRange{3,1}(3), arrayLUTRange{3,2}(3), numStepsInLUT(3))');
        arrayColLUT{4} = arrayLUTRange{4,1};
        
        %move through each LUT
        iOutputLUT = 1;
        for iLUT = 1:numLUTs,

            numSteps = numStepsInLUT(iLUT);

            %identify variable pairs with a correlation within the specified
            %range
            for iStep = 1:numSteps,

                %apply the LUT
                numMaxVal = arrayDataToColLUT{iLUT}(iStep);
                numMinVal = arrayDataToColLUT{iLUT}(iStep+1);
                if iLUT < 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData <= numMaxVal) & (arrayInputData > numMinVal));
                elseif iLUT == 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData <= numMaxVal) & (arrayInputData >= numMinVal));
                end

                %and colour the output array
                for iPixel = 1:length(arrayPixelRow),
                    arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayColLUT{iLUT}(iStep,:);
                end

                %output the LUT to be displayed next to the heat map
                arrayLUTcoloured(iOutputLUT,1,:) = arrayColLUT{iLUT}(iStep,:);
                iOutputLUT = iOutputLUT+1;

            end

        end

        %specify NaN values to another colour
        arrayIsNaN = isnan(arrayInputData);
        [arrayPixelRow, arrayPixelCol] = find(arrayIsNaN);
        for iPixel = 1:length(arrayPixelRow),
            arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayNaNLUT;
        end
        
        %specify change points on the LUT for display
        structLUTInfo(1).LUTPix = 0.5;
        structLUTInfo(1).DispValue = num2str(max(arrayInputData(:)), '%0.3f');
        structLUTInfo(2).LUTPix = numStepsInLUT(1)+0.5;
        structLUTInfo(2).DispValue = num2str(structDataInfo.HighThresh, '%0.3f');
        structLUTInfo(3).LUTPix = sum(numStepsInLUT(1:3))+0.5;
        structLUTInfo(3).DispValue = num2str(structDataInfo.LowThresh, '%0.3f');
        structLUTInfo(4).LUTPix = sum(numStepsInLUT(1:4))+0.5;
        structLUTInfo(4).DispValue = '0';
        
    elseif strcmp(structDataInfo.Type,'mutinfo') || strcmp(structDataInfo.Type,'data'),
        
        %calculate the offset around the 'zero' range
        numZeroOffSet = structDataInfo.Thresh/(numStepsInLUT(2)+0.5);

        %specify the correlation values for the colour LUT
        arrayMutInfoToColLUT = cell(numLUTs,1);
        arrayMutInfoToColLUT{1} = linspace(max(arrayInputData(:)),structDataInfo.Thresh, numStepsInLUT(1)+1);
        arrayMutInfoToColLUT{2} = linspace(structDataInfo.Thresh, numZeroOffSet, numStepsInLUT(2)+1);
        arrayMutInfoToColLUT{3} = [numZeroOffSet 0];
        
        %specify the colours within the LUT
        arrayColLUT = cell(numLUTs,1);
        arrayColLUT{1} = cat(2, linspace(arrayLUTRange{1,1}(1), arrayLUTRange{1,2}(1), numStepsInLUT(1))', linspace(arrayLUTRange{1,1}(2), arrayLUTRange{1,2}(2), numStepsInLUT(1))', linspace(arrayLUTRange{1,1}(3), arrayLUTRange{1,2}(3), numStepsInLUT(1))');
        arrayColLUT{2} = cat(2, linspace(arrayLUTRange{2,1}(1), arrayLUTRange{2,2}(1), numStepsInLUT(2))', linspace(arrayLUTRange{2,1}(2), arrayLUTRange{2,2}(2), numStepsInLUT(2))', linspace(arrayLUTRange{2,1}(3), arrayLUTRange{2,2}(3), numStepsInLUT(2))');
        arrayColLUT{3} = arrayLUTRange{3,1};
                
        %move through each LUT
        iOutputLUT = 1;
        for iLUT = 1:numLUTs,

            numSteps = numStepsInLUT(iLUT);

            %identify variable pairs with a correlation within the specified
            %range
            for iStep = 1:numSteps,

                %apply the LUT
                numMaxVal = arrayMutInfoToColLUT{iLUT}(iStep);
                numMinVal = arrayMutInfoToColLUT{iLUT}(iStep+1);
                if iLUT < 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData <= numMaxVal) & (arrayInputData > numMinVal));
                elseif iLUT == 3,
                    [arrayPixelRow, arrayPixelCol] = find((arrayInputData <= numMaxVal) & (arrayInputData >= numMinVal));
                end

                %and colour the output array
                for iPixel = 1:length(arrayPixelRow),
                    arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayColLUT{iLUT}(iStep,:);
                end

                %output the LUT to be displayed next to the heat map
                arrayLUTcoloured(iOutputLUT,1,:) = arrayColLUT{iLUT}(iStep,:);
                iOutputLUT = iOutputLUT+1;

            end

        end

        %specify NaN values to another colour
        arrayIsNaN = isnan(arrayInputData);
        [arrayPixelRow, arrayPixelCol] = find(arrayIsNaN);
        for iPixel = 1:length(arrayPixelRow),
            arrayHeatMapcoloured(arrayPixelRow(iPixel), arrayPixelCol(iPixel), :) = arrayNaNLUT;
        end
        
        %specify change points on the LUT for display
        structLUTInfo(1).LUTPix = 0.5;
        structLUTInfo(1).DispValue = num2str(max(arrayInputData(:)), '%0.3f');
        structLUTInfo(2).LUTPix = numStepsInLUT(1)+0.5;
        structLUTInfo(2).DispValue = num2str(structDataInfo.Thresh, '%0.3f');
        structLUTInfo(3).LUTPix = sum(numStepsInLUT(1:3))+0.5;
        structLUTInfo(3).DispValue = '0';
        
    end


end

