import os
import csv
import numpy as np
import matplotlib.pyplot as plt

__author__ = 'jcursons'
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is used to analyse melanoma (SKCM) data from The Cancer Genome Atlas (TCGA) project, in particular to
#  visualise the expression of micro-RNAs and putative target mRNAs across clinical samples. This analysis produces
#  a subset of the results presented in:
#   MC Andrews/J Cursons, DG Hurley, M Anaka, JS Cebon, A Behren, EJ Crampin (2016). Systems analysis identifies
#    miR-29b regulation of invasiveness in melanoma. BMC Molecular Cancer, (Accepted Nov 2016).
#   http://dx.doi.org/doi-to-be-assigned
#
# Further details on this project can be found on the GitHub repository:
#  http://github.com/uomsystemsbiology/LMMEL-miR-miner
#
# For further information, please contact:
#  Dr. Miles Andrews
#   MD Anderson Cancer Centre
#   miles.andrews (at) onjcri.org.au
#
#  Dr. Joe Cursons
#   Bioinformatics Division, Walter and Eliza Hall Institute of Medical Research, Australia
#   cursons.j (at) wehi.edu.au
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script contains three functions which are used to extract micro-RNA (miR) and mRNA abundances from the TCGA SKCM
#  data, and match data between patients, defined within:
# class AnalyseTCGA
#  --> function extractMicRNAData: combine miR data across multiple files and export a dictionary containing arrays and
#                                   pointers
#  --> function extractMessRNAData: combine mRNA data across multiple files and export a dictionary containing arrays
#                                   and pointers
#  --> function matchMicAndMessRNAData: examine TCGA sample barcodes between the miR and mRNA data, and create an
#                                       array/vector which points between related samples.
#
# These functions are executed by the appended script which plots a number of specified associations.
#  * hsa-mir-29b-1 and hsa-mir-29b-2 abundance are plotted against the abundance of LAMC1, PPIC and LASP1
#       --> Fig4_miR29b_TCGAplots.png/Fig4_miR29b_TCGAplots.eps
#  * the abundance for a number of miRs are plotted against the abundance of mRNAs identified as putative targets
#     through the analysis of the LM-MEL cell line panel
#       --> Fig2_TCGAPlots.png
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class AnalyseTCGA:

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # extractMicRNAData
    #
    # Inputs:
    #   strInMicRNADataDir: a string containing the file path for the TCGA SKCM micro-RNA data
    #   flagPerformExtraction: a Boolean flag which specifies whether the data need to be extracted (and saved). This
    #                           parameter must be set to True for the initial execution, but can then be set to False to
    #                           improve run time on subsequent execution.
    # Outputs:
    #   A dictionary containing label vectors and the data array within matched key/value pairs:
    #       'data': a 2D np.float64 array containing the miR transcript abundance, with micro-RNAs in rows and TCGA
    #               samples in columns
    #       'miRLabels': a list containing the names of micro-RNAs, matched to the data array
    #       'obsLabels': a list containing the names of TCGA samples, matched to the data array
    # # # # # # # # # #
    def extractMicRNAData(strInMicRNADataDir, flagPerformExtraction):

        # specify a name for the processed miR data file (so that data extraction only needs to be performed once, to
        #  improve runtime on subsequent executions)
        strOutputSaveFile = 'processedMicRNAData'
        # identify the base directory path for the miR data by performing an rsplit of the string around the path
        #  separator
        # TODO: write this using os.path
        strTCGABaseDir = strInMicRNADataDir.rsplit("\\", 1)[0]

        # if flagPerformExtraction=False but the processed data file is not present, overwrite this setting and provide
        #  a warning to the user
        if np.bitwise_and((not os.path.exists(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))),
                          (not flagPerformExtraction)):
            print('warning: the pre-processed miR data cannot be found at ' +
                  os.path.exists(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz'))) +
                  ' and flagPerformExtraction=False; \n' +
                  '          forcing flagPerformExtraction=True, this will increase the run time')
            flagPerformExtraction=True

        # if specified, perform the miR data extraction
        if flagPerformExtraction:
            print('Attempting to extract miR data from ' + strInMicRNADataDir)
            # extract the name of all files within the specified folder
            arrayInputFileNames = os.listdir(strInMicRNADataDir)
            numFiles = len(arrayInputFileNames)

            # step through every file and perform a pre-index of the miRs
            listObsFileLabels = []
            print('... Performing initial extraction of miR names')
            for iFile in range(numFiles):

                # process the string for the observation file label
                strExpectedSuffix = '.isoform.quantification.txt'
                numExpectedSuffixLength = len(strExpectedSuffix)
                strMicRNAFileName = arrayInputFileNames[iFile]
                numFileNameStringLength = len(strMicRNAFileName)

                # determine the observation file label from the file name
                if strMicRNAFileName[(numFileNameStringLength-numExpectedSuffixLength):] == strExpectedSuffix:
                    strObsFileLabel = strMicRNAFileName[0:(numFileNameStringLength-numExpectedSuffixLength)]
                    listObsFileLabels.append(strObsFileLabel)
                else:
                    print('The miR data filename was not in the expected format to extract the label')

                # load the file
                filePointer = open(os.path.join(strInMicRNADataDir, arrayInputFileNames[iFile]))
                arrayFile = csv.reader(filePointer, delimiter='\t')

                # convert to a list and determine its length
                arrayFileAsList = list(arrayFile)
                numRows = len(arrayFileAsList)

                # extract all mature miR labels; exclude the header row (start at index 1)
                listTempMatureMicRNAList = []
                for iRow in range(1,numRows):
                    arrayRow = arrayFileAsList[iRow]
                    # the miR processing status/form is in the 5th column
                    if arrayRow[5][0:6] == 'mature':
                        listTempMatureMicRNAList.append(arrayRow[0])

                # convert all miRs within this file to a (unique) set
                setTempMatureMicRNAs = set(listTempMatureMicRNAList)

                # close the csv
                filePointer.close()

                if iFile == 0:
                    # assign the current set to be the full set
                    setOutputMicRNAs = setTempMatureMicRNAs
                else:
                    # combine the sets to define a more comprehensive set
                    setOutputMicRNAs = setTempMatureMicRNAs | setOutputMicRNAs

            # extract back to a list and determine the length
            listMicRNANames = list(setOutputMicRNAs)
            numMicRNAs = len(listMicRNANames)

            print('... Extracting miR names and data across all files')

            # create an array which stores the combined data
            arrayCombinedMicRNAData = np.zeros((numMicRNAs, numFiles), dtype=np.float64)

            # it can take a bit of time to extract the miR data, so create a vector for progress reporting
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
                        listTempMatureMicRNAList.append(arrayRow[0])

                        # the miR name is in the first column
                        strMicRNA = arrayRow[0]
                        # determine the output index
                        numOutputMicRNAIndex = listMicRNANames.index(strMicRNA)

                        # the RPKM-normalised count data are in the fourth column
                        numDataObs = np.float64(arrayRow[3])

                        # sum over all mature counts
                        arrayCombinedMicRNAData[numOutputMicRNAIndex, iFile] = \
                            arrayCombinedMicRNAData[numOutputMicRNAIndex, iFile] + numDataObs

            # save the output arrays using numpy (as this is quicker than reading csv files back in)
            np.savez(os.path.join(strTCGABaseDir, strOutputSaveFile),
                     arrayMicRNANames=listMicRNANames,
                     arrayCombinedMicRNAData=arrayCombinedMicRNAData,
                     arrayObsFileLabels=listObsFileLabels
                     )

        else:
            if os.path.exists(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz'))):
                print('Loading the processed miR data from ' +
                      os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))

                # load the data from the specified files
                npzfile = np.load(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))

                # extract the arrays/vectors from the dict for output
                arrayCombinedMicRNAData = npzfile['arrayCombinedMicRNAData']
                listMicRNANames = npzfile['arrayMicRNANames']
                listObsFileLabels = npzfile['arrayObsFileLabels']

            else:
                # in theory the function shouldn't get here because I check the presence of the data and the value for
                #  flagPerformExtraction above..
                print('Cannot load the miR data, ' + os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')) +
                      ' does not exist, change flagPerformExtraction')

        return {'data': arrayCombinedMicRNAData, 'miRLabels': listMicRNANames, 'obsLabels': listObsFileLabels}

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # extractMessRNAData
    #
    # Inputs:
    #   strInMessRNADataDir: a string containing the file path for the TCGA SKCM mRNA data
    #   flagPerformExtraction: a Boolean flag which specifies whether the data need to be extracted (and saved). This
    #                           parameter must be set to True for the initial execution, but can then be set to False to
    #                           improve run time on subsequent execution.
    # Outputs:
    #   A dictionary containing label vectors and the data array within matched key/value pairs:
    #       'data': a 2D np.float64 array containing the miR transcript abundance, with mRNAs in rows and TCGA
    #               samples in columns
    #       'mRNALabels': a list containing the HGNC symbol for mRNAs, matched to the data array
    #       'obsLabels': a list containing the names of TCGA samples, matched to the data array
    # # # # # # # # # #
    def extractMessRNAData(strInMessRNADataDir, flagPerformExtraction):

        strOutputSaveFile = 'processedMessRNAData'
        #TODO: write this using os.path
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

                # extract the mRNA data filename
                strMessRNAFileName = arrayInputFileNames[iFile]
                numFileNameStringLength = len(strMessRNAFileName)

                # check that the filename prefix and suffix are in the format expected
                flagIsExpectedPrefix = (strMessRNAFileName[:(numExpectedPrefixLength)] == strExpectedPrefix)
                flagIsExpectedSuffix = (strMessRNAFileName[(numFileNameStringLength-numExpectedSuffixLength):] == strExpectedSuffix)
                if flagIsExpectedPrefix and flagIsExpectedSuffix:
                    # extract the observation file label and output
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

                # move through the first file and extract all of the miR names into a list
                listMessRNANames = []
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

                    # extract the mRNA name from before the pipe/bar/"|"
                    stringMessRNAName = stringNameData[0:numPipePos]
                    stringMessRNANum = stringNameData[(numPipePos+1):]

                    # and output
                    listMessRNANames.append(stringMessRNAName)
                    arrayMessRNAEntrezIDs[iRow-1] = int(stringMessRNANum)

                # create an array which stores the combined data
                arrayCombinedMessRNAData = np.zeros((numMessRNAs, numFiles), dtype=float)

                # it can take a bit of time to extract the miR data, so create a vector for progress reporting
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

                    # load the text file
                    with open(os.path.join(strInMessRNADataDir, arrayInputFileNames[iFile])) as filePointer:
                        # extract the file contents using csv.reader
                        arrayFile = csv.reader(filePointer, delimiter='\t')
                        # convert the file contents to a list
                        arrayFileAsList = list(arrayFile)

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
            np.savez(os.path.join(strTCGABaseDir, strOutputSaveFile), arrayMessRNANames=listMessRNANames, arrayCombinedMessRNAData=arrayCombinedMessRNAData, arrayObsFileLabels=arrayObsFileLabels)


        else:
            if os.path.exists(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz'))):
                # load the data from the specified files
                print('Loading the processed mRNA data from ' + os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))

                npzfile = np.load(os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')))

                listMessRNANames = npzfile['arrayMessRNANames']
                arrayCombinedMessRNAData = npzfile['arrayCombinedMessRNAData']
                arrayObsFileLabels = npzfile['arrayObsFileLabels']

            else:
                print('Cannot load the mRNA data, ' + os.path.join(strTCGABaseDir, (strOutputSaveFile + '.npz')) + ' does not exist, change flagPerformExtraction')

        return {'data': arrayCombinedMessRNAData, 'mRNALabels': listMessRNANames, 'obsLabels': arrayObsFileLabels}



    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # matchMicAndMessRNAData
    #
    # Inputs:
    #   strInTCGAFolderPath: a string containing the file path for the TCGA SKCM mRNA data
    #   arrayInMicRNABarcodes: an array of miR sample barcodes
    #   arrayInMessRNAUUIDs: an array of mRNA sample UUIDs
    #
    # Outputs:
    #   arrayForMicRNACorrMessRNAIndex: a vector which lists the miR vector index which matches the mRNA sample
    #                                    at the same position in the barcode vector
    # # # # # # # # # #
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Execute the code with specified variables. Here we are attempting to generate two sets of figures which are integrated
#  with other results to produce figures within the published manuscript.
#
# The first section of code plots the abundance of hsa-mir-29b-1 and hsa-mir-29b-2 abundance against the abundance of
#  putative targets identified by the MATLAB script, in particular: LAMC1, PPIC and LASP1
#       --> Fig4_miR29b_TCGAplots.png/Fig4_miR29b_TCGAplots.eps
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# specify the micro-RNAs within the TCGA SKCM data; note that these data map to immature micro-RNAs (mir-x rather than
#  miR-x), and thus there are two version of hsa-mir-29b (-1 and -2)
# miRBase accession IDs corresponding to these targets are:
#  	MI0000105 = hsa-mir-29b-1
#  	MI0000107 = hsa-mir-29b-2
listMicRNAsToPlot = ['hsa-mir-29b-1', 'hsa-mir-29b-2']
numMicRNAsToPlot = len(listMicRNAsToPlot)

# specify HGNC symbols for the mRNAs within the TCGA SKCM data
listMessRNAsToPlot = ['LAMC1', 'PPIC', 'LASP1']
numMessRNAsToPlot = len(listMessRNAsToPlot)

# control data processing (i.e. load intermediate processed files to decrease run time)
flagPerformMicRNAExtraction = False
flagPerformMessRNAExtraction = False

# define the file system location of the input files
strTCGASKCMBasePath = 'C:\\db\\tcga\\skcm'
strMicRNADataPath = os.path.join(strTCGASKCMBasePath, 'miR')
strMessRNADataPath = os.path.join(strTCGASKCMBasePath, 'mRNA')

# define the file system location of the output files; just output to the location of this script
strOutputFolder = os.path.dirname(os.path.realpath(__file__))
# check that the folder exists, if not, create
if not os.path.exists(strOutputFolder):
    os.makedirs(strOutputFolder)

# extract the miR and mRNA abundance data, together with associated pointer arrays
structMicRNAData = AnalyseTCGA.extractMicRNAData(strMicRNADataPath, flagPerformMicRNAExtraction)
structMessRNAData = AnalyseTCGA.extractMessRNAData(strMessRNADataPath, flagPerformMessRNAExtraction)

# use the pointer arrays and TCGA sample map file to combine miR and mRNA abundance data
arrayMicMessRNAPairedIndices = AnalyseTCGA.matchMicAndMessRNAData(strTCGASKCMBasePath,
                                                                  structMicRNAData['obsLabels'],
                                                                  structMessRNAData['obsLabels'])
numObs = len(arrayMicMessRNAPairedIndices)

# specify variables for creating the output figure
numAnnotationFontSize = 10
numMaxYTicks = 3
numMaxXTicks = 3

# create a multi-panel figure for the final output
handleFig, arrayAxesHandles = plt.subplots(2,3)
handleFig.set_size_inches(w=6.3, h=4.3)

# step through each specified miR
for iMicRNA in range(numMicRNAsToPlot):

    # match the string to the pointer array
    stringMicRNAToTest = listMicRNAsToPlot[iMicRNA]
    # and extract the data index (while checking the number of matches)
    if stringMicRNAToTest in structMicRNAData['miRLabels']:
        numMicRNAIndex = np.where(structMicRNAData['miRLabels'] == stringMicRNAToTest)
        if len(numMicRNAIndex) == 1:
            # only a single transcript has been matched
            numMicRNAIndex = numMicRNAIndex[0]
        else:
            # multiple transcripts have been matched
            numMicRNAIndex = numMicRNAIndex[0]
            print('warning: multiple indices exist for ' + stringMicRNAToTest + ', using the first index')
    else:
        # no transcript has been matched
        numMicRNAIndex = 0
        print('warning: ' + stringMicRNAToTest + ' cannot be found within the TCGA miR list')

    # step through each mRNA specified
    for iMessRNA in range(numMessRNAsToPlot):

        # match the mRNA string to the pointer array
        stringMessRNAToTest = listMessRNAsToPlot[iMessRNA]
        # and extract the data index (while checking the number of matches)
        if stringMessRNAToTest in structMessRNAData['mRNALabels']:
            numMessRNAIndex = np.where(structMessRNAData['mRNALabels'] == stringMessRNAToTest)
            if len(numMessRNAIndex) == 1:
                # only a single transcript has been matched
                numMessRNAIndex = numMessRNAIndex[0]
            else:
                # multiple transcripts have been matched
                numMessRNAIndex = numMessRNAIndex[0]
                print('warning: multiple indices exist for ' + stringMessRNAToTest + ', using the first index')
        else:
            # no transcript has been matched
            numMessRNAIndex = 0
            print('warning: ' + stringMessRNAToTest + ' cannot be found within the TCGA mRNA list')

        # check that the data have been found within the data matrices
        if (numMessRNAIndex > 0) and (numMicRNAIndex > 0):
            # create matched vectors for paired miR:mRNA data/observations
            arrayMicRNAVector = np.zeros(numObs, dtype=float)
            arrayMessRNAVector = np.zeros(numObs, dtype=float)

            # step through each patient/observation
            for iObs in range(len(arrayMicMessRNAPairedIndices)):
                # extract appropriate data points into the matched vectors
                arrayMicRNAVector[iObs] = structMicRNAData['data'][numMicRNAIndex, iObs]
                arrayMessRNAVector[iObs] = structMessRNAData['data'][numMessRNAIndex,arrayMicMessRNAPairedIndices[iObs]]

            # from the paired observations, calculate the Pearson's correlation
            numTCGALogDataCorr = np.corrcoef(np.log2(arrayMicRNAVector), np.log2(arrayMessRNAVector))[0, 1]
            # note that corrcoef calculates the correlation between all pairwise combinations of inputs, including
            #  self-correlations (== 1) along the diagonal; thus extract the [0,1] data point

            # create a scatter plot displaying the matched data
            arrayAxesHandles[iMicRNA,iMessRNA].scatter(np.log2(arrayMicRNAVector), np.log2(arrayMessRNAVector),
                                                       s=5, c='0.8', edgecolors='0.3', alpha=0.6)

            # label the x- and y-axes
            arrayAxesHandles[iMicRNA,iMessRNA].set_xlabel(stringMicRNAToTest + '\nabundance',
                                                          fontsize=numAnnotationFontSize)
            arrayAxesHandles[iMicRNA,iMessRNA].set_ylabel(stringMessRNAToTest + ' abundance',
                                                          fontsize=numAnnotationFontSize)

            # calculate a linear best fit to the data
            structFit = np.polyfit(np.log2(arrayMicRNAVector), np.log2(arrayMessRNAVector), 1)
            funcFit = np.poly1d(structFit)

            # overlay the line of best fit
            numMinX = min(np.log2(arrayMicRNAVector))
            numMaxX = max(np.log2(arrayMicRNAVector))
            arrayAxesHandles[iMicRNA,iMessRNA].plot([numMinX, numMaxX], funcFit([numMinX, numMaxX]), 'r-', lw=2)

            # tidy up the x-axis tick locations
            arrayXTickLoc = plt.MaxNLocator(numMaxXTicks)
            arrayAxesHandles[iMicRNA,iMessRNA].xaxis.set_major_locator(arrayXTickLoc)

            # tidy up the y-axis tick locations
            arrayYTickLoc = plt.MaxNLocator(numMaxYTicks)
            arrayAxesHandles[iMicRNA,iMessRNA].yaxis.set_major_locator(arrayYTickLoc)

            # specify the tick font size
            arrayAxesHandles[iMicRNA,iMessRNA].tick_params(labelsize=numAnnotationFontSize)

            # hide the right and top spines
            arrayAxesHandles[iMicRNA,iMessRNA].spines['right'].set_visible(False)
            arrayAxesHandles[iMicRNA,iMessRNA].spines['top'].set_visible(False)

            # only show ticks on the left and bottom spines
            arrayAxesHandles[iMicRNA,iMessRNA].yaxis.set_ticks_position('left')
            arrayAxesHandles[iMicRNA,iMessRNA].xaxis.set_ticks_position('bottom')

            # extract the x- and y-axis data limits
            [numMinXAxis, numMaxXAxis] = arrayAxesHandles[iMicRNA,iMessRNA].get_xlim()
            [numMinYAxis, numMaxYAxis] = arrayAxesHandles[iMicRNA,iMessRNA].get_ylim()

            # use these to calculate an appropriate position for overlaying the calculated correlation
            arrayAxesHandles[iMicRNA,iMessRNA].text((numMinXAxis + (numMaxXAxis-numMinXAxis)*0.8),
                                                    (numMinYAxis + (numMaxYAxis-numMinYAxis)*0.95),
                                                    ('$r_P$ = ' + "{0:.3f}".format(numTCGALogDataCorr)),
                                                    color='r', horizontalalignment='center',
                                                    verticalalignment='center', fontsize=numAnnotationFontSize)

# use the tight_layout function to minimise non-functional white space
handleFig.tight_layout()

# output the figure as a png file
handleFig.savefig(os.path.join(strOutputFolder, 'Fig4_miR29b_TCGAplots.png'), format='png', dpi=300)

# output the figure as an eps file
handleFig.savefig(os.path.join(strOutputFolder, 'Fig4_miR29b_TCGAplots.eps'), format='eps', dpi=600)

# close the output figure
plt.close(handleFig)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# The second section of code plots the abundance of various micro-RNAs against the abundance of predicted/putative
#  target mRNAs (identified by the MATLAB script). For specific details please refer to 'listRelsToPlot'.
#  putative targets identified by the MATLAB script, in particular: LAMC1, PPIC and LASP1
#       --> Fig2_TCGAPlots.png
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# create a list of lists; where embedded lists contain [miR, mRNA] pairs for putative relationships
listRelsToPlot = [['hsa-let-7b', 'LIN28B'],
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
                  ['hsa-mir-500a', 'AXL'],
                  ['hsa-mir-30b', 'FOXD1'],
                  ['hsa-mir-185', 'NRP1'],
                  ['hsa-mir-146a', 'TCF4'],
                  ['hsa-mir-211', 'TCF4'],
                  ['hsa-mir-125b-1*', 'HPS4'],
                  ['hsa-let-7e', 'CAPN3'],
                  ['hsa-mir-125b-1', 'GYG2'],
                  ['hsa-mir-125b-2', 'GYG2'],
                  ['hsa-mir-199b', 'MBP'],
                  ['hsa-mir-502', 'SOX9'],
                  ['hsa-mir-222', 'SOX10']]
numRelsToPlot = len(listRelsToPlot)

# specify plotting parameters for individual plots
numAnnotationFontSize = 8
numMaxYTicks = 3
numMaxXTicks = 3

# specify parameters for creating the subplot array
numPlotsOverX = 6
numPlotsOverY = 8

# create a multi-panel figure for the final output
handleFig, arrayAxesHandles = plt.subplots(numPlotsOverY, numPlotsOverX)
handleFig.set_size_inches(w=66/numPlotsOverY, h=60/numPlotsOverX)

# step through each relationship/miR:mRNA pair
for iRel in range(numRelsToPlot):

    # calculate the local subplot array position
    numPlotRow = np.int(np.floor(iRel/numPlotsOverX)*2)
    numPlotCol = np.int(iRel - ((numPlotRow/2)*numPlotsOverX))

    # extract the miR string
    stringMicRNAToTest = listRelsToPlot[iRel][0]
    if stringMicRNAToTest in structMicRNAData['miRLabels']:
        # match to the full list of TCGA SKCM miRs and extract the index, while checking the number of matches
        numMicRNAIndex = np.where(structMicRNAData['miRLabels'] == stringMicRNAToTest)
        if len(numMicRNAIndex) == 1:
            # only a single match
            numMicRNAIndex = numMicRNAIndex[0]
        else:
            # multiple matches identified
            numMicRNAIndex = numMicRNAIndex[0]
            print('warning: multiple indices exist for ' + stringMicRNAToTest + ', using the first index')
    else:
        # no match identified
        numMicRNAIndex = 0
        print('warning: ' + stringMicRNAToTest + ' cannot be found within the TCGA miR list')

    # extract the mRNA string
    stringMessRNAToTest = listRelsToPlot[iRel][1]
    if stringMessRNAToTest in structMessRNAData['mRNALabels']:
        # match to the full list of TCGA SKCM mRNAs and extract the index, while checking the number of matches
        numMessRNAIndex = np.where(structMessRNAData['mRNALabels'] == stringMessRNAToTest)
        if len(numMessRNAIndex) == 1:
            # only a single match
            numMessRNAIndex = numMessRNAIndex[0]
        else:
            # multiple matches identified
            numMessRNAIndex = numMessRNAIndex[0]
            print('warning: multiple indices exist for ' + stringMessRNAToTest + ', using the first index')
    else:
        # no match identified
        numMessRNAIndex = 0
        print('warning: ' + stringMessRNAToTest + ' cannot be found within the TCGA mRNA list')

    # if data are found for both the miR and mRNA specified
    if (numMessRNAIndex > 0) and (numMicRNAIndex > 0):
        # create matched vectors for the data
        arrayMicRNAVector = np.zeros(numObs, dtype=float)
        arrayMessRNAVector = np.zeros(numObs, dtype=float)

        # step through each patient/observation
        for iObs in range(len(arrayMicMessRNAPairedIndices)):
            # extract the appropriate data into the matched vectors
            arrayMicRNAVector[iObs] = structMicRNAData['data'][numMicRNAIndex, iObs]
            arrayMessRNAVector[iObs] = structMessRNAData['data'][numMessRNAIndex,arrayMicMessRNAPairedIndices[iObs]]

        # identify any samples with non-zero abundance for both specified transcripts
        arrayPairedObs = (arrayMessRNAVector > 0) & (arrayMicRNAVector > 0)

        # from the paired observations, calculate the correlations; as noted above, extract the number in [0,1] from
        #  the corrcoef output array
        numTCGALogDataCorr = np.corrcoef(np.log2(arrayMicRNAVector[arrayPairedObs]),
                                         np.log2(arrayMessRNAVector[arrayPairedObs]))[0, 1]

        # prune the hsa- prefix and match back
        if stringMicRNAToTest[0:7] == 'hsa-let':
            stringMicRNAToDisp = stringMicRNAToTest[4:]
        elif stringMicRNAToTest[0:7] == 'hsa-mir':
            stringMicRNAToDisp = 'miR-' + stringMicRNAToTest[8:]

        # create a scatter plot displaying the matched data
        arrayAxesHandles[numPlotRow, numPlotCol].scatter(np.log2(arrayMicRNAVector[arrayPairedObs]),
                                                         np.log2(arrayMessRNAVector[arrayPairedObs]),
                                                         s=5, c='0.8', edgecolors='0.3', alpha=0.7)

        # specify the x- and y- axis labels
        arrayAxesHandles[numPlotRow, numPlotCol].set_xlabel(stringMicRNAToDisp,
                                                            fontsize=numAnnotationFontSize)
        arrayAxesHandles[numPlotRow, numPlotCol].set_ylabel(stringMessRNAToTest,
                                                            fontsize=numAnnotationFontSize)

        # calculate a linear best fit to the data
        structFit = np.polyfit(np.log2(arrayMicRNAVector[arrayPairedObs]),
                               np.log2(arrayMessRNAVector[arrayPairedObs]), 1)
        funcFit = np.poly1d(structFit)

        # overlay the line of best fit
        numMinX = min(np.log2(arrayMicRNAVector[arrayPairedObs]))
        numMaxX = max(np.log2(arrayMicRNAVector[arrayPairedObs]))
        arrayAxesHandles[numPlotRow, numPlotCol].plot([numMinX, numMaxX], funcFit([numMinX, numMaxX]), 'r-')

        # tidy up the x-axis tick locations
        arrayXTickLoc = plt.MaxNLocator(numMaxXTicks)
        arrayAxesHandles[numPlotRow, numPlotCol].xaxis.set_major_locator(arrayXTickLoc)

        # tidy up the y-axis tick locations
        arrayYTickLoc = plt.MaxNLocator(numMaxYTicks)
        arrayAxesHandles[numPlotRow, numPlotCol].yaxis.set_major_locator(arrayYTickLoc)

        # specify the tick label font size
        arrayAxesHandles[numPlotRow, numPlotCol].tick_params(labelsize=numAnnotationFontSize)

        # extract the x- and y-axis limits
        [numMinXAxis, numMaxXAxis] = arrayAxesHandles[numPlotRow, numPlotCol].get_xlim()
        [numMinYAxis, numMaxYAxis] = arrayAxesHandles[numPlotRow, numPlotCol].get_ylim()

        # use the axis limits to specify the position of the calculated Pearson's correlation
        arrayAxesHandles[numPlotRow, numPlotCol].text((numMinXAxis + (numMaxXAxis-numMinXAxis)*0.8),
                                                      (numMinYAxis + (numMaxYAxis-numMinYAxis)*0.9),
                                                      ('$r_P$ = ' + "{0:.2f}".format(numTCGALogDataCorr)),
                                                      color='r', horizontalalignment='center',
                                                      verticalalignment='center', fontsize=numAnnotationFontSize)

        # hide the axis for the 'plot above' (spacing to combine the TCGA and LM-MEL plots)
        arrayAxesHandles[(numPlotRow-1), numPlotCol].axis('off')

        # hide the right and top spines
        arrayAxesHandles[numPlotRow, numPlotCol].spines['right'].set_visible(False)
        arrayAxesHandles[numPlotRow, numPlotCol].spines['top'].set_visible(False)

        # only show ticks on the left and bottom spines
        arrayAxesHandles[numPlotRow, numPlotCol].yaxis.set_ticks_position('left')
        arrayAxesHandles[numPlotRow, numPlotCol].xaxis.set_ticks_position('bottom')

# use the tight_layout function to minimise non-functional white space
handleFig.tight_layout()

# output the figure as a png file
handleFig.savefig(os.path.join(strOutputFolder, 'Fig2_TCGAPlots'), ext='png', close=True, dpi=300)

# close the output figure
plt.close(handleFig)
