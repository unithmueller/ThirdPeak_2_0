#Manage Python invocation from Matlab
import sys
import os
sys.path.append(os.path.normpath(os.getcwd()+"\\Scripts\\PredictData"))
from TemporalSegandFPMatlab import TemporalSegandFPMatlab
from PredictChangepointByPretrainedNN import PredictDataUsingPretrainedChangepointLSM
from PredictProteinbyMLP import PredictProteinbyMLP

#%% DifferentAnalyisModes
def performTempSeqAndFingerprinting(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew):
    CleanFingerprintData, results3, referenceList = TemporalSegandFPMatlab(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, False)
    #print("return value of size" + str(len(results3)) + "and type" + str(type(results3)))
    return results3

def performTempSeqFingerprintingChangepointEst(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, DownstreamModelSavePath):
    CleanFingerprintData, results3, referenceList = TemporalSegandFPMatlab(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, False)
    test_outputs_list, test_probs_list, test_changepoint_pred = PredictDataUsingPretrainedChangepointLSM(CleanFingerprintData, DownstreamModelSavePath, filePath)
    return results3, test_outputs_list, test_probs_list, test_changepoint_pred

def performTempSeqFingerprintingClassification(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, ClassificModel, ConfThresh):
    CleanFingerprintData, results3, referenceList = TemporalSegandFPMatlab(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, False)
    y_pred = PredictProteinbyMLP(ClassificModel, CleanFingerprintData, ConfThresh, filePath)
    return results3, referenceList, y_pred

def performTempSeqOnly(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew):
    ensemblScore, ensemblePred = TemporalSegandFPMatlab(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, True)
    return ensemblScore, ensemblePred



#%% main managing function
def manageDeepSPTCallsfromMatlab(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, DownstreamModelSavePath, ConfThresh, modeselection):
    print("")
    #Controls what to do with data from matlab, calls the underlying python functions and sends the data back to matlab in the end
    #modeselection = 0
    if modeselection == 0:
        print("mode 0")
        fingerprintedTracksforMatlab = performTempSeqAndFingerprinting(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew)
        print("got file to return:" + str(type(fingerprintedTracksforMatlab)) + "with length" + str(len(fingerprintedTracksforMatlab)))   
        return fingerprintedTracksforMatlab, [], [], []
    elif modeselection == 1:
        print("mode 1")
        results3, test_outputs_list, test_probs_list, test_changepoint_pred = performTempSeqFingerprintingChangepointEst(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, DownstreamModelSavePath)
        print("got files to return:" + str(type(results3)) + "and " + str(type(test_outputs_list)) + "and"  + str(type(test_probs_list))) 
        return results3, test_outputs_list, test_probs_list, test_changepoint_pred
    elif modeselection == 2:
        print("mode 2")
        results3, referenceList, y_pred = performTempSeqFingerprintingClassification(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, DownstreamModelSavePath, ConfThresh)
        print("got files to return:" + str(type(results3)) + "and " + str(type(referenceList)) + "and"  + str(type(y_pred))) 
        return results3, referenceList, y_pred
    elif modeselection == 3:
        print("mode 3")
        ensemblScore, ensemblePred = performTempSeqOnly(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew)
        print("got files to return:" + str(type(ensemblScore)) + "and " + str(type(ensemblePred))) 
        return ensemblScore, ensemblePred
          
#%% Run the file

results, test_outputs_list, test_probs_list, test_changepoint_pred = manageDeepSPTCallsfromMatlab(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew, DownstreamModelSavePath, ConfThresh, modeselection)