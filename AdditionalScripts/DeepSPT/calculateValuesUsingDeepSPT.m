function calculatedTracks = calculateValuesUsingDeepSPT(Trackfiles, filePath, timestep, dimension, windowSize)
    %Calculate the temporal segmentation and Fingerprints of the given
    %tracks. Saves files to their original filepaths to retrieve them
    %later for further processing if necessary.
    
    %% Import options for the usage of certain models during the processing
    functionName = 'calculateValuesUsingDeepSPT.m';
    settingsTextPath = which(functionName);
    settingsTextPath = settingsTextPath(1:end-(size(functionName,2)));
    temporalSegmentation3D_NN = join([settingsTextPath "temporalSegmentationNN3D_used.txt"],"");
    temporalSegmentation2D_NN = join([settingsTextPath "temporalSegmentationNN2D_used.txt"],"");
    hMM_used = join([settingsTextPath "hMM_used.txt"],"");
    LSM_used = join([settingsTextPath "Downstream_LSM_used.txt"],"");
    MLP_used = join([settingsTextPath "Dowstream_MLP_used.txt"],"");
    mode_used = join([settingsTextPath "mode.txt"],"");
    
    %Read the files
    temporalSegmentation3D_NN = readlines(temporalSegmentation3D_NN);
    temporalSegmentation2D_NN = readlines(temporalSegmentation2D_NN);
    hMM_used = readlines(hMM_used);
    LSM_used = readlines(LSM_used);
    MLP_used = readlines(MLP_used);
    mode_used = readlines(mode_used);


    %put the paths into a list as it make some problems with the
    %delimiter...
    temporalSegmentation3D_NN = split(temporalSegmentation3D_NN,["/", "\"]);
    temporalSegmentation2D_NN = split(temporalSegmentation2D_NN,["/", "\"]);
    hMM_used = split(hMM_used,["/", "\"]);
    LSM_used = split(LSM_used,["/", "\"]);
    MLP_used = split(MLP_used,["/", "\"]);

    temporalSegmentation3D_NN = join(temporalSegmentation3D_NN,"/");
    temporalSegmentation2D_NN = join(temporalSegmentation2D_NN,"/");
    hMM_used = join(hMM_used,"/");
    LSM_used = join(LSM_used,"/");
    MLP_used = join(MLP_used,"/");
    
    if dimension == 2
        usedTmpSegNN = temporalSegmentation2D_NN;
    else
        usedTmpSegNN = temporalSegmentation3D_NN;
    end
    
    usedHMM = hMM_used;
    generateNew = 0;
    ConfThresh = 0.6;
    modeselection = double(mode_used);

    if modeselection == 0
        DownstreamModelSavePath = "";
    elseif modeselection == 1
        DownstreamModelSavePath = LSM_used;
    elseif modeselction == 2
        DownstreamModelSavePath = MLP_used;
    else
        DownstreamModelSavePath = "";
    end

    %% Prepare tracks 
    numberOfFiles = size(Trackfiles,1);
    FileLengths = [];
    calculatedTracks = cell(numberOfFiles,2);
    convertedFiles = py.list();
    convertedPaths = py.list();

    for i = 1:numberOfFiles
        tmpLen = size(Trackfiles{i,1},1);
        FileLengths(i) = tmpLen;
        tmp = py.numpy.array(Trackfiles{i,1});
        tmpName = Trackfiles{i,2};
        calculatedTracks{i,2} = tmpName;
        convertedFiles.append(tmp);
        convertedPaths.append(tmpName);
    end

    %% Send to python
    processedTracks = pyrunfile("ManageDeepSPTCallsfromMatlab.py", "results", matlabTracks=convertedFiles, fileNames=convertedPaths, filePath=filePath, dimension=dimension, timestep=timestep, windowSize=windowSize, usedTmpSegNN=usedTmpSegNN, usedHMM=usedHMM, generateNew=generateNew, DownstreamModelSavePath=DownstreamModelSavePath, ConfThresh=ConfThresh, modeselection=modeselection);
    
    %% Repack the data for matlab
    
    tmpResults = cell(processedTracks);
    tmpResults = cellfun(@double, tmpResults, 'UniformOutput', false);
    tmpResults = cell2mat(tmpResults);
    tmpResults = reshape(tmpResults, [], 43);
    
    idxCounter = 1;

    for i = 1:numberOfFiles
        tmp = Trackfiles{i,1};
        tmpLen = FileLengths(i);
        tmp(:,23:65) = tmpResults(idxCounter:(idxCounter+tmpLen-1),:);
        calculatedTracks{i,1} = tmp;
        idxCounter = idxCounter+tmpLen;
    end
end