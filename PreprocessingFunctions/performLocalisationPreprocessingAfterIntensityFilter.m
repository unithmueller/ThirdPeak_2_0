function driftCorrLocs = performLocalisationPreprocessingAfterIntensityFilter(SaveFolderPath, SaveFolderName, loadedIntensityLocs, BeadLocations, options, saveCSV)
    %Main function for the preprocessing of the localisation data.
   %Input: SaveFolderPath: Path on the Disk to save to
   %SaveFolderName: Name of the subfolder generated in the SaveFolderPath
   %to save to
   %ImportSettingsStruct: Contains information of the file types being used
   %FileLocations: List of paths to the files that will be processed
   %Options: Options from the GUI for filtering
   %Output: Filtered Localisations
    %% set save location and go there
    cd(SaveFolderPath);
    counter = 1;
    foldername = iterateSaveFoldername(SaveFolderName, counter);
    cd(foldername);
    
    %% set variables for saving
    %intensityLocs = {};

    %% sideload intensityfilterd Locs
    intensityLocs = loadedIntensityLocs;
   
    %% perform the drift correction
    if options.drift.Performdrift == 1
        %do drift correction
        if options.drift.ReferenceAvailable == 1
            %calculate by drift of a bead
            driftCorrLocs = performPreprocessingDriftCorrectionWithBead(BeadLocations, intensityLocs);
        else
            %calculate by the data
            driftCorrLocs = performPreprocessingDriftCorrectionwoBead(intensityLocs);
        end
    else
        driftCorrLocs = intensityLocs;
    end
    %% save the data to disk
    %save the localisations in mat and for swift
    datatosave = intensityLocs;
    save("intensityFilteredLocalisations.mat","datatosave");
     if saveCSV
        mkdir intensityFilteredCSV
        tmpPath = pwd + "/intensityFilteredCSV";
        for f = 1:size(datatosave,1)
            tmpFileData = datatosave{f,1};
            [~, tmpFileName, ~] = fileparts(datatosave{f,2});
            writematrix(tmpFileData, fullfile(tmpPath, tmpFileName, ".csv"));
        end
    end
    datatosave = driftCorrLocs;
    save("driftCorrectedLocalisations.mat","datatosave");
    if saveCSV
        mkdir driftCorrectedLocalisationsCSV
        tmpPath = pwd + "/driftCorrectedLocalisationsCSV";
        for f = 1:size(datatosave,1)
            tmpFileData = datatosave{f,1};
            [~, tmpFileName, ~] = fileparts(datatosave{f,2});
            writematrix(tmpFileData, fullfile(tmpPath, tmpFileName, ".csv"));
        end
    end
    %save the filter settings
    setvalues = options;
    save("PropertiesUsed.mat","setvalues");
    writestruct(setvalues,"PropertiesUsed.xml");
    %save for swift
    SavePeaksAsSwift(pwd, driftCorrLocs, size(intensityLocs,1), 1);
end