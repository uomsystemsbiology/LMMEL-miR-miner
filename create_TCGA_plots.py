import os
import csv
import numpy as np
import matplotlib.pyplot as plt

__author__ = 'jcursons'

class AnalyseTCGA:

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # function which moves through the specified miR data folder and extracts the required information
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def extractMicRNAData(strInMicRNADataDir, flagPerformExtraction):

        strOutputSaveFile = 'processedMicRNAData'
        strTCGABaseDir = strInMicRNADataDir.rsplit("\\", 1)[0]

        if flagPerformExtraction:
            print('Attempting to extract miR data from ' + strInMicRNADataDir)
            # extract the name of all files within the specified folder
            arrayInputFileNames = os.listdir(strInMicRNADataDir)
            numFiles = len(arrayInputFileNames)

            # step through every file and perform a pre-index of the miRs
            arrayObsFileLabels = []
            print('... Performing initial extraction of miR names')
            for iFile in range(numFiles):

                # process the string for the observation file label
                strExpectedSuffix = '.isoform.quantification.txt'
                numExpectedSuffixLength = len(strExpectedSuffix)
                strMicRNAFileName = arrayInputFileNames[iFile]
                numFileNameStringLength = len(strMicRNAFileName)

                if strMicRNAFileName[(numFileNameStringLength-numExpectedSuffixLength):] == strExpectedSuffix:
                    strObsFileLabel = strMicRNAFileName[0:(numFileNameStringLength-numExpectedSuffixLength)]
                    arrayObsFileLabels.append(strObsFileLabel)
                else:
                    print('The miR data filename was not in the expected format to extract the label')

                # load the file
                filePointer = open(os.path.join(strInMicRNADataDir, arrayInputFileNames[iFile]))
                arrayFile = csv.reader(filePointer, delimiter='\t')

                # convert to a list and determine its length
                arrayFileAsList = list(arrayFile)
                numRows = len(arrayFileAsList)

                # extract all mature miR labels; exclude the header row (start at index 1)
                arrayTempMatureMicRNAList = []
                for iRow in range(1,numRows):
                    arrayRow = arrayFileAsList[iRow]
                    # the miR processing status/form is in the 5th column
                    if arrayRow[5][0:6] == 'mature':
                        arrayTempMatureMicRNAList.append(arrayRow[0])

                # convert all miRs within this file to a (unique) set
                arrayTempMatureMicRNASet = set(arrayTempMatureMicRNAList)

                # close the csv
                filePointer.close()

                if iFile == 0:
                    # assign the current set to be the full set
                    arrayOutputMicRNASet = arrayTempMatureMicRNASet
                else:
                    # combine the sets to define a more comprehensive set
                    arrayOutputMicRNASet = arrayTempMatureMicRNASet | arrayOutputMicRNASet

            # extract back to a list and determine the length
            arrayMicRNANames = list(arrayOutputMicRNASet)
            numMicRNAs = len(arrayMicRNANames)

            print('... Extracting miR names and data across all files')

            # create an array which stores the combined data
            arrayCombinedMicRNAData = np.zeros((numMicRNAs, numFiles), dtype=np.float64)

            numProgDivisions = 20
            numProgCounter = numFiles/numProgDivisions
            numProgOut = 100/numProgDivisions

            # step through all files and extract the data
            for iFile in range(numFiles):

                # output script progress for some of the larger data sets
                if iFile > numProgCounter:
                    print('...... ' + str(numProgOut) + '% through ' + str(numFiles) + ' files')
                    numProgCounter = numProgCounter + numFiles/numProgDivisions
                    numProgOut = numProgOut + 100/numProgDivisions

                # open the file
                filePointer = open(os.path.join(strInMicRNADataDir, arrayInputFileNames[iFile]))
                # extract the tab-delimited contents into an array using csv
                arrayFile = csv.reader(filePointer, delimiter='\t')
                # convert the file contents to a list for subsequent indexing
                arrayFileAsList = list(arrayFile)
                numRows = len(arrayFileAsList)
                #and close the file
                filePointer.close()

                # move through the file, skipping the header row ([0])
                for iRow in range(1, numRows):

                    # extract the row
                    arrayRow = arrayFileAsList[iRow]

                    # the miR processing status/form is in the 5th column
                    if arrayRow[5][0:6] == 'mature':
                        arrayTempMatureMicRNAList.append(arrayRow[0])

                        # the miR name is in the first column
                        strMicRNA = arrayRow[0]
                        # determine the output index
                        numOutputMicRNAIndex = arrayMicRNANames.index(strMicRNA)

                        # the RPKM-normalised count data are in the fourth column
                        numDataObs = np.float64(arrayRow[3])

                        # sum over all mature counts
                        arrayCombinedMicRNAData[numOutputMicRNAIndex, iFile] = arrayCombinedMicRNAData[numOutputMicRNAIndex, iFile] + numDataObs

            # save the output arrays using numpy (as this is quicker than reading csv files back in)
            np.savez(os.path.join(strTCGABaseDir, strOutputSaveFile), arrayMicRNANames=arrayMicRNANames, arrayCombinedMicRNAData=arrayCombinedMicRNAData, arrayObsFileLabels=arrayObsFileLabels)

        else:
            if os.path.exists(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz'))):
                # load the data from the specified files
                print('Loading the processed miR data from ' + os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))

                npzfile = np.load(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))

                arrayCombinedMicRNAData = npzfile['arrayCombinedMicRNAData']
                arrayMicRNANames = npzfile['arrayMicRNANames']
                arrayObsFileLabels = npzfile['arrayObsFileLabels']

            else:
                print('Cannot load the miR data, ' + os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')) + ' does not exist, change flagPerformExtraction')

        return {'data': arrayCombinedMicRNAData, 'miRLabels': arrayMicRNANames, 'obsLabels': arrayObsFileLabels}

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # function which moves through the specified mRNA data folder and extracts the required information
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def extractMessRNAData(strInMessRNADataDir, flagPerformExtraction):

        strOutputSaveFile = 'processedMessRNAData'
        strTCGABaseDir = strInMessRNADataDir.rsplit("\\", 1)[0]

        if flagPerformExtraction:
            print('Attempting to extract mRNA data from ' + strInMessRNADataDir)
            # extract the name of all files within the specified folder
            arrayInputFileNames = os.listdir(strInMessRNADataDir)
            numFiles = len(arrayInputFileNames)

            # step through every file and check that they all contain the same number of entries/rows
            arrayNumMessRNAs = []
            arrayObsFileLabels = []

            print('... Performing initial extraction of file names and checking file sizes')
            for iFile in range(numFiles):

                # process the string for the observation file label
                strExpectedSuffix = '.rsem.genes.normalized_results'
                strExpectedPrefix = 'unc.edu.'
                numExpectedSuffixLength = len(strExpectedSuffix)
                numExpectedPrefixLength = len(strExpectedPrefix)

                strMessRNAFileName = arrayInputFileNames[iFile]
                numFileNameStringLength = len(strMessRNAFileName)

                flagIsExpectedPrefix = (strMessRNAFileName[:(numExpectedPrefixLength)] == strExpectedPrefix)
                flagIsExpectedSuffix = (strMessRNAFileName[(numFileNameStringLength-numExpectedSuffixLength):] == strExpectedSuffix)

                if flagIsExpectedPrefix and flagIsExpectedSuffix:
                    strObsFileLabel = strMessRNAFileName[(numExpectedPrefixLength):(numFileNameStringLength-numExpectedSuffixLength)]
                    arrayObsFileLabels.append(strObsFileLabel)
                else:
                    print('The mRNA data filename was not in the expected format to extract the label')

                # load the file and determine its length
                filePointer = open(os.path.join(strInMessRNADataDir, arrayInputFileNames[iFile]))
                arrayFile = csv.reader(filePointer, delimiter='\t')
                numRows = len(list(arrayFile))

                # remove the header line to get the number of miRs
                arrayNumMessRNAs.append(numRows-1)
                filePointer.close()

            # convert to a set to produce a unique list
            setNumMessRNAs = set(arrayNumMessRNAs)

            if len(setNumMessRNAs) == 1:
                print('... Extracting mRNA names and data across all files')
                # the number of miRs across all files is equal, so use the first value; subtract one for the header row
                numMessRNAs = arrayNumMessRNAs[0]

                # move through the first file and extract all of the miR names into an array
                arrayMessRNANames = []
                arrayMessRNAEntrezIDs = np.zeros((numMessRNAs, 1), dtype=np.uint32)

                filePointer = open(os.path.join(strInMessRNADataDir, arrayInputFileNames[0]))
                arrayFile = csv.reader(filePointer, delimiter='\t')

                # convert to a list and step through every row
                arrayFileAsList = list(arrayFile)
                filePointer.close()

                for iRow in range(1, numMessRNAs):
                    arrayRow = arrayFileAsList[iRow]
                    # extract the mRNA name and Entrez ID from the first column of each row
                    stringNameData = arrayRow[0]
                    numPipePos = stringNameData.find('|')

                    stringMessRNAName = stringNameData[0:numPipePos]
                    stringMessRNANum = stringNameData[(numPipePos+1):]

                    arrayMessRNANames.append(stringMessRNAName)
                    arrayMessRNAEntrezIDs[iRow-1] = int(stringMessRNANum)

                # create an array which stores the combined data
                arrayCombinedMessRNAData = np.zeros((numMessRNAs, numFiles), dtype=float)

                numProgDivisions = 20
                numProgCounter = numFiles/numProgDivisions
                numProgOut = 100/numProgDivisions

                # step through all files and extract the data
                for iFile in range(numFiles):

                    # output script progress for some of the larger data sets
                    if iFile > numProgCounter:
                        print('...... ' + str(numProgOut) + '% through ' + str(numFiles) + ' files')
                        numProgCounter = numProgCounter + numFiles/numProgDivisions
                        numProgOut = numProgOut + 100/numProgDivisions

                    filePointer = open(os.path.join(strInMessRNADataDir, arrayInputFileNames[iFile]))
                    arrayFile = csv.reader(filePointer, delimiter='\t')

                    arrayFileAsList = list(arrayFile)
                    filePointer.close()

                    for iRow in range(1, numMessRNAs):
                        # the first element of the row is the mRNA name/number, and the second element is the normalised count data
                        arrayRow = arrayFileAsList[iRow]
                        # extract the mRNA name and Entrez ID from the first column of each row
                        stringNameData = arrayRow[0]
                        numPipePos = stringNameData.find('|')

                        numMessRNANum = int(stringNameData[(numPipePos+1):])

                        numDataObs = float(arrayRow[1])

                        # check that the correct miR is being extracted
                        if numMessRNANum == arrayMessRNAEntrezIDs[iRow-1]:
                            arrayCombinedMessRNAData[iRow-1, iFile] = numDataObs

                        else:
                            print("The order of the mRNAs is not as expected, the data matching code needs to be generalised")

            else:
                # spit out a warning
                print("The file " " does not contain the expected number of mRNAs")

            # save the output arrays using numpy (as this is quicker than reading csv files back in)
            np.savez(os.path.join(strTCGABaseDir, strOutputSaveFile), arrayMessRNANames=arrayMessRNANames, arrayCombinedMessRNAData=arrayCombinedMessRNAData, arrayObsFileLabels=arrayObsFileLabels)


        else:
            if os.path.exists(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz'))):
                # load the data from the specified files
                print('Loading the processed mRNA data from ' + os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))

                npzfile = np.load(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))

                arrayMessRNANames = npzfile['arrayMessRNANames']
                arrayCombinedMessRNAData = npzfile['arrayCombinedMessRNAData']
                arrayObsFileLabels = npzfile['arrayObsFileLabels']

            else:
                print('Cannot load the mRNA data, ' + os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')) + ' does not exist, change flagPerformExtraction')

        return {'data': arrayCombinedMessRNAData, 'mRNALabels': arrayMessRNANames, 'obsLabels': arrayObsFileLabels}



    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # function which matches miR and mRNA data files
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def matchMicAndMessRNAData(strInTCGAFolderPath, arrayInMicRNABarcodes, arrayInMessRNAUUIDs):
        # load the file
        filePointer = open(os.path.join(strInTCGAFolderPath, 'FILE_SAMPLE_MAP.txt'))
        arrayFile = csv.reader(filePointer, delimiter='\t')

        arrayFileAsList = list(arrayFile)

        # create an array which contains the TSS and participant sections of the TCGA barcode from the sample map file,
        #  and a corresponding array which contains file names
        arrayTCGABarcode = []
        arrayCorrFile = []

        # create an output array
        arrayForMicRNACorrMessRNAIndex = np.zeros((len(arrayInMicRNABarcodes), 1), dtype='uint16')

        # remove the last three fields of the TCGA barcode
        arrayProcMicRNABarcodes = []
        for stringInMicRNABarcode in arrayInMicRNABarcodes:
            arrayProcMicRNABarcodes.append(stringInMicRNABarcode.rsplit('-', 2)[0])

        for iRow in range(1,len(arrayFileAsList)-1):
            stringRow = arrayFileAsList[iRow]

            # the TCGA barcode is located within the second column
            stringBarcode = stringRow[1]
            arrayTCGABarcode.append(stringBarcode)

            # the file name is located within the first column
            stringFileName = stringRow[0]
            arrayCorrFile.append(stringFileName)


        # step through all of the MessRNA UUIDs and attempt to match back to a barcode
        for iFile in range(len(arrayInMessRNAUUIDs)):
            stringMessRNAUUID = arrayInMessRNAUUIDs[iFile]
            stringFileToFind = 'unc.edu.' + stringMessRNAUUID + '.rsem.genes.normalized_results'

            # check if this file is listed in the 'corresponding filenames' array
            if stringFileToFind in arrayCorrFile:
                numMessRNAUUIDIndexMatch = arrayCorrFile.index(stringFileToFind)

                # if so, extract the TCGA barcode for the mRNA file
                stringMessRNABarcode = arrayTCGABarcode[numMessRNAUUIDIndexMatch]

                #remove the final 3 fields also
                stringMessRNABarcodeToMatch = stringMessRNABarcode.rsplit('-',2)[0]

                # and attempt to match this to a MicRNA barcodes being input
                if stringMessRNABarcodeToMatch in arrayProcMicRNABarcodes:
                    numMicRNAIndexMatch = arrayProcMicRNABarcodes.index(stringMessRNABarcodeToMatch)

                    if type(numMicRNAIndexMatch) is not list:
                        arrayForMicRNACorrMessRNAIndex[numMicRNAIndexMatch] = iFile
                    else:
                        print('warning: multiple mRNA files have been matched to the miR file through TCGA barcode')
                        arrayForMicRNACorrMessRNAIndex[numMicRNAIndexMatch] = iFile

                else:
                    print('warning: the mRNA barcode cannot be matched to a MicRNA barcode')

            else:
                print('warning: cannot find the mRNA file name in the sample map lists')

        return arrayForMicRNACorrMessRNAIndex





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# execute the code with specified variables
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

arrayMicRNAsToPlot = ['hsa-mir-29b-1', 'hsa-mir-29b-2']

numMicRNAsToPlot = len(arrayMicRNAsToPlot)

arrayMessRNAsToPlot = ['LAMC1', 'PPIC', 'LASP1']

numMessRNAsToPlot = len(arrayMessRNAsToPlot)

# control data processing (i.e. load intermediate processed files to decrease run time)
flagPerformMicRNAExtraction = False
flagPerformMessRNAExtraction = False

# define the file system location of the input files
strTCGASKCMBasePath = 'C:\\db\\tcga\\skcm'
strMicRNADataPath = strTCGASKCMBasePath + '\\miR'
strMessRNADataPath = strTCGASKCMBasePath + '\\mRNA'

# define the file system location of the output files
strOutputFolder = os.path.dirname(os.path.realpath(__file__))
# check that the folder exists, if not, create
if not os.path.exists(strOutputFolder):
    os.makedirs(strOutputFolder)

# extract the miR and mRNA abundance data, together with associated pointer arrays
structMicRNAData = AnalyseTCGA.extractMicRNAData(strMicRNADataPath, flagPerformMicRNAExtraction)
structMessRNAData = AnalyseTCGA.extractMessRNAData(strMessRNADataPath, flagPerformMessRNAExtraction)

# use the pointer arrays and TCGA sample map file to combine miR and mRNA abundance data
arrayMicMessRNAPairedIndices = AnalyseTCGA.matchMicAndMessRNAData(strTCGASKCMBasePath, structMicRNAData['obsLabels'], structMessRNAData['obsLabels'])
numObs = len(arrayMicMessRNAPairedIndices)

numAnnotationFontSize = 8
numMaxYTicks = 3
numMaxXTicks = 4

# create a multi-panel figure for the final output
handleFig, arrayAxesHandles = plt.subplots(2,3)
handleFig.set_size_inches(6.3, 4.3)

for iMicRNA in range(numMicRNAsToPlot):

    # match the string to the pointer array
    stringMicRNAToTest = arrayMicRNAsToPlot[iMicRNA]
    if stringMicRNAToTest in structMicRNAData['miRLabels']:
        numMicRNAIndex = np.where(structMicRNAData['miRLabels'] == stringMicRNAToTest)
        if len(numMicRNAIndex) == 1:
            numMicRNAIndex = numMicRNAIndex[0]
        else:
            numMicRNAIndex = numMicRNAIndex[0]
            print('warning: multiple indices exist for ' + stringMicRNAToTest + ', using the first index')
    else:
        numMicRNAIndex = 0
        print('warning: ' + stringMicRNAToTest + ' cannot be found within the TCGA miR list')

    # step through each miR-mRNA pair
    for iMessRNA in range(numMessRNAsToPlot):

        # match the mRNA string to the pointer array
        stringMessRNAToTest = arrayMessRNAsToPlot[iMessRNA]
        if stringMessRNAToTest in structMessRNAData['mRNALabels']:
            numMessRNAIndex = np.where(structMessRNAData['mRNALabels'] == stringMessRNAToTest)
            if len(numMessRNAIndex) == 1:
                numMessRNAIndex = numMessRNAIndex[0]
            else:
                numMessRNAIndex = numMessRNAIndex[0]
                print('warning: multiple indices exist for ' + stringMessRNAToTest + ', using the first index')
        else:
            numMessRNAIndex = 0
            print('warning: ' + stringMessRNAToTest + ' cannot be found within the TCGA mRNA list')

        # create arrays with paired observations
        if (numMessRNAIndex > 0) and (numMicRNAIndex > 0):
            arrayMicRNAVector = np.zeros(numObs, dtype=float)
            arrayMessRNAVector = np.zeros(numObs, dtype=float)

            for iObs in range(len(arrayMicMessRNAPairedIndices)):

                arrayMicRNAVector[iObs] = structMicRNAData['data'][numMicRNAIndex, iObs]
                arrayMessRNAVector[iObs] = structMessRNAData['data'][numMessRNAIndex,arrayMicMessRNAPairedIndices[iObs]]

            # from the paired observations, calculate the correlations
            numTCGALogDataCorr = np.corrcoef(np.log2(arrayMicRNAVector), np.log2(arrayMessRNAVector))[0, 1]

            # and output appropriate plots
            arrayAxesHandles[iMicRNA,iMessRNA].scatter(np.log2(arrayMicRNAVector), np.log2(arrayMessRNAVector), s=5, c='0.8', edgecolors='0.3')
            arrayAxesHandles[iMicRNA,iMessRNA].set_xlabel(stringMicRNAToTest + ' abundance', fontsize=numAnnotationFontSize)
            arrayAxesHandles[iMicRNA,iMessRNA].set_ylabel(stringMessRNAToTest + ' abundance', fontsize=numAnnotationFontSize)

            # calculate a linear best fit to the data
            structFit = np.polyfit(np.log2(arrayMicRNAVector), np.log2(arrayMessRNAVector), 1)
            funcFit = np.poly1d(structFit)

            # overlay the line of best fit
            numMinX = min(np.log2(arrayMicRNAVector))
            numMaxX = max(np.log2(arrayMicRNAVector))
            arrayAxesHandles[iMicRNA,iMessRNA].plot([numMinX, numMaxX], funcFit([numMinX, numMaxX]), 'r-')

            # tidy up the x-axis tick locations
            arrayXTickLoc = plt.MaxNLocator(numMaxXTicks)
            arrayAxesHandles[iMicRNA,iMessRNA].xaxis.set_major_locator(arrayXTickLoc)

            arrayYTickLoc = plt.MaxNLocator(numMaxYTicks)
            arrayAxesHandles[iMicRNA,iMessRNA].yaxis.set_major_locator(arrayYTickLoc)

            arrayAxesHandles[iMicRNA,iMessRNA].tick_params(labelsize=numAnnotationFontSize)

            [numMinXAxis, numMaxXAxis] = arrayAxesHandles[iMicRNA,iMessRNA].get_xlim()
            [numMinYAxis, numMaxYAxis] = arrayAxesHandles[iMicRNA,iMessRNA].get_ylim()

            # hide the right and top spines
            arrayAxesHandles[iMicRNA,iMessRNA].spines['right'].set_visible(False)
            arrayAxesHandles[iMicRNA,iMessRNA].spines['top'].set_visible(False)

            # only show ticks on the left and bottom spines
            arrayAxesHandles[iMicRNA,iMessRNA].yaxis.set_ticks_position('left')
            arrayAxesHandles[iMicRNA,iMessRNA].xaxis.set_ticks_position('bottom')

            arrayAxesHandles[iMicRNA,iMessRNA].text((numMinXAxis + (numMaxXAxis-numMinXAxis)*0.8), (numMinYAxis + (numMaxYAxis-numMinYAxis)*0.95),
                                            ('$r_P$ = ' + "{0:.3f}".format(numTCGALogDataCorr)), color='r',
                                            horizontalalignment='center', verticalalignment='center', fontsize=numAnnotationFontSize)

handleFig.tight_layout()

handleFig.savefig(os.path.join(strOutputFolder, 'FigX_miR29b_TCGAplots'), ext='png', close=True, dpi=300)
plt.close(handleFig)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# execute the code with specified variables
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

arrayRelsToPlot = [['hsa-let-7b', 'LIN28B'],
                   ['hsa-mir-30b', 'RUNX2'],
                   ['hsa-mir-30b', 'SERPINE1'],
                   ['hsa-mir-125b-1', 'IRF4'],
                   ['hsa-mir-125b-2', 'IRF4'],
                   ['hsa-mir-211', 'TGFBR2'],
                   ['hsa-mir-29b-1', 'CDK6'],
                   ['hsa-mir-29b-2', 'CDK6'],
                   ['hsa-mir-29b-1', 'COL4A1'],
                   ['hsa-mir-29b-2', 'COL4A1'],
                   ['hsa-mir-29b-1', 'PDGFC'],
                   ['hsa-mir-29b-2', 'PDGFC'],
                   ['hsa-mir-17', 'CYBRD1'],
                   ['hsa-mir-30d', 'CPE'],
                   ['hsa-mir-30d', 'JUN'],
                   ['hsa-mir-98', 'HIC2'],
                   ['hsa-mir-146a', 'NRAS'],
                   ['hsa-mir-185', 'NRP1'],
                   ['hsa-mir-185', 'SRPX2'],
                   ['hsa-mir-199b', 'MBP'],
                   ['hsa-mir-211', 'ANXA1'],
                   ['hsa-mir-211', 'TCF4'],
                   ['hsa-mir-222', 'CHKA'],
                   ['hsa-mir-222', 'SOX10']]

# associations that don't really hold within the TCGA data:
#                    ['hsa-let-7e', 'CHD7'],
#                    ['hsa-mir-19b-1', 'RASA1'],
#                    ['hsa-mir-19b-2', 'RASA1'],
#                    ['hsa-mir-21', 'SEMA6A'],
#                    ['hsa-mir-23b', 'BCL2'],
#                    ['hsa-mir-29a', 'C5orf15'],
#                    ['hsa-mir-29a', 'COL4A1'],
#                    ['hsa-mir-29b-1', 'CDC42'],
#                    ['hsa-mir-29b-2', 'CDC42'],
#                    ['hsa-mir-30d', 'GNAI1'],
#                    ['hsa-mir-34a', 'NOTCH1'],
#                    ['hsa-mir-34a', 'FOSL1'],
#                    ['hsa-mir-93', 'BMPR2'],
#                    ['hsa-mir-106b', 'BMPR2'],
#                    ['hsa-mir-125b-1', 'ZEB2'],
#                    ['hsa-mir-125b-2', 'ZEB2'],
#                    ['hsa-mir-125b-2', 'IRF4'],
#                    ['hsa-mir-148a', 'SERPINE1'],
#                    ['hsa-mir-148a', 'NRP1'],
#                    ['hsa-mir-182', 'EDNRB'],
#                    ['hsa-mir-221', 'ZFYVE16'],
#                    ['hsa-mir-221', 'RRAGD'],
#                    ['hsa-mir-222', 'BAMBI'],
#                    ['hsa-mir-222', 'TGFBR3'],
#                    ['hsa-mir-224', 'MBP'],
#                    ['hsa-mir-340', 'THBS1'],
#                    ['hsa-mir-340', 'PLAU'],
#                    ['hsa-mir-340', 'BMP2'],
#                    ['hsa-mir-361', 'RAB35'],
#                    ['hsa-mir-500a', 'THBS1'],
#                    ['hsa-mir-532', 'THBS1']

numRelsToPlot = len(arrayRelsToPlot)

numAnnotationFontSize = 8

numMaxYTicks = 3
numMaxXTicks = 3

numPlotsOverX = 6
numPlotsOverY = 8

# create a multi-panel figure for the final output
handleFig, arrayAxesHandles = plt.subplots(numPlotsOverY, numPlotsOverX)
handleFig.set_size_inches(66/numPlotsOverY, 60/numPlotsOverX)

for iRel in range(numRelsToPlot):

    numPlotRow = np.floor(iRel/numPlotsOverX)*2
    numPlotCol = iRel - ((numPlotRow/2)*numPlotsOverX)

    # match the string to the pointer array
    stringMicRNAToTest = arrayRelsToPlot[iRel][0]
    if stringMicRNAToTest in structMicRNAData['miRLabels']:
        numMicRNAIndex = np.where(structMicRNAData['miRLabels'] == stringMicRNAToTest)
        if len(numMicRNAIndex) == 1:
            numMicRNAIndex = numMicRNAIndex[0]
        else:
            numMicRNAIndex = numMicRNAIndex[0]
            print('warning: multiple indices exist for ' + stringMicRNAToTest + ', using the first index')
    else:
        numMicRNAIndex = 0
        print('warning: ' + stringMicRNAToTest + ' cannot be found within the TCGA miR list')

    # match the mRNA string to the pointer array
    stringMessRNAToTest = arrayRelsToPlot[iRel][1]
    if stringMessRNAToTest in structMessRNAData['mRNALabels']:
        numMessRNAIndex = np.where(structMessRNAData['mRNALabels'] == stringMessRNAToTest)
        if len(numMessRNAIndex) == 1:
            numMessRNAIndex = numMessRNAIndex[0]
        else:
            numMessRNAIndex = numMessRNAIndex[0]
            print('warning: multiple indices exist for ' + stringMessRNAToTest + ', using the first index')
    else:
        numMessRNAIndex = 0
        print('warning: ' + stringMessRNAToTest + ' cannot be found within the TCGA mRNA list')

    # create arrays with paired observations
    if (numMessRNAIndex > 0) and (numMicRNAIndex > 0):
        arrayMicRNAVector = np.zeros(numObs, dtype=float)
        arrayMessRNAVector = np.zeros(numObs, dtype=float)

        for iObs in range(len(arrayMicMessRNAPairedIndices)):

            arrayMicRNAVector[iObs] = structMicRNAData['data'][numMicRNAIndex, iObs]
            arrayMessRNAVector[iObs] = structMessRNAData['data'][numMessRNAIndex,arrayMicMessRNAPairedIndices[iObs]]

        if any(arrayMessRNAVector == 0):
            arrayMessRNAVector[arrayMessRNAVector == 0] = 0.04

        if any(arrayMicRNAVector == 0):
            arrayMicRNAVector[arrayMicRNAVector == 0] = 0.04

        arrayPairedObs = (arrayMessRNAVector > 0.04) & (arrayMicRNAVector > 0.04)

        # from the paired observations, calculate the correlations
        numTCGALogDataCorr = np.corrcoef(np.log2(arrayMicRNAVector[arrayPairedObs]), np.log2(arrayMessRNAVector[arrayPairedObs]))[0, 1]

        if stringMicRNAToTest[0:7] == 'hsa-let':
            stringMicRNAToDisp = stringMicRNAToTest[4:]
        elif stringMicRNAToTest[0:7] == 'hsa-mir':
            stringMicRNAToDisp = 'miR-' + stringMicRNAToTest[8:]

        # and output appropriate plots
        arrayAxesHandles[numPlotRow, numPlotCol].scatter(np.log2(arrayMicRNAVector), np.log2(arrayMessRNAVector), s=5, c='0.8', edgecolors='0.3')
        arrayAxesHandles[numPlotRow, numPlotCol].set_xlabel(stringMicRNAToDisp, fontsize=numAnnotationFontSize)
        arrayAxesHandles[numPlotRow, numPlotCol].set_ylabel(stringMessRNAToTest, fontsize=numAnnotationFontSize)

        # calculate a linear best fit to the data
        structFit = np.polyfit(np.log2(arrayMicRNAVector[arrayPairedObs]), np.log2(arrayMessRNAVector[arrayPairedObs]), 1)
        funcFit = np.poly1d(structFit)

        # overlay the line of best fit
        numMinX = min(np.log2(arrayMicRNAVector))
        numMaxX = max(np.log2(arrayMicRNAVector))
        arrayAxesHandles[numPlotRow, numPlotCol].plot([numMinX, numMaxX], funcFit([numMinX, numMaxX]), 'r-')

        # tidy up the x-axis tick locations
        arrayXTickLoc = plt.MaxNLocator(numMaxXTicks)
        arrayAxesHandles[numPlotRow, numPlotCol].xaxis.set_major_locator(arrayXTickLoc)

        arrayYTickLoc = plt.MaxNLocator(numMaxYTicks)
        arrayAxesHandles[numPlotRow, numPlotCol].yaxis.set_major_locator(arrayYTickLoc)

        arrayAxesHandles[numPlotRow, numPlotCol].tick_params(labelsize=numAnnotationFontSize)

        [numMinXAxis, numMaxXAxis] = arrayAxesHandles[numPlotRow, numPlotCol].get_xlim()
        [numMinYAxis, numMaxYAxis] = arrayAxesHandles[numPlotRow, numPlotCol].get_ylim()

        arrayAxesHandles[numPlotRow, numPlotCol].text((numMinXAxis + (numMaxXAxis-numMinXAxis)*0.8), (numMinYAxis + (numMaxYAxis-numMinYAxis)*0.9),
                                        ('$r_P$ = ' + "{0:.3f}".format(numTCGALogDataCorr)), color='r',
                                        horizontalalignment='center', verticalalignment='center', fontsize=numAnnotationFontSize)

        # hide the axis for the 'plot above' (spacing to combine the TCGA and LM-MEL plots)
        arrayAxesHandles[(numPlotRow-1), numPlotCol].axis('off')

        # hide the right and top spines
        arrayAxesHandles[numPlotRow, numPlotCol].spines['right'].set_visible(False)
        arrayAxesHandles[numPlotRow, numPlotCol].spines['top'].set_visible(False)

        # only show ticks on the left and bottom spines
        arrayAxesHandles[numPlotRow, numPlotCol].yaxis.set_ticks_position('left')
        arrayAxesHandles[numPlotRow, numPlotCol].xaxis.set_ticks_position('bottom')

handleFig.tight_layout()

handleFig.savefig(os.path.join(strOutputFolder, 'Fig2_TCGAPlots'), ext='png', close=True, dpi=300)
plt.close(handleFig)
