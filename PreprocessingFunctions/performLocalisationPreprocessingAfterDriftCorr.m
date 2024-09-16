function driftCorrLocs = performLocalisationPreprocessingAfterDriftCorr(SaveFolderPath, SaveFolderName, loadedDriftLocs, options, saveCSV)
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

    %%Load the data
   driftCorrLocs = loadedDriftLocs;

    %% save the data to disk
    %save the localisations in mat and for swift
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
    SavePeaksAsSwift(pwd, driftCorrLocs, size(driftCorrLocs,1), 1);
end