function driftCorrLocs = performLocalisationPreprocessing(SaveFolderPath, SaveFolderName, ImportSettingsStruct, FileLocations, BeadLocations, options, saveCSV)
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
    %% grab file list and load files depending on settings
    LocalisationData = loadLocswithFileList(ImportSettingsStruct, FileLocations);
    try %optional
        BeadLocations = loadLocswithFileList(ImportSettingsStruct, BeadLocations);
    catch
    end

    %make sure there are no negative values for XY
    for i = 1:size(LocalisationData,1)
        tmpdata = LocalisationData{i,1};
        xdat = tmpdata(:,2);
        ydat = tmpdata(:,3);
        xdat(xdat<0) = 0;
        ydat(ydat<0) = 0;
        tmpdata(:,2:3) = [xdat,ydat];
        LocalisationData{i,1} = tmpdata;
    end


    clear xdat ydat tmpdata;
    %% set variables for saving
    maskedLocs = cell(size(LocalisationData,1)*20,2); %allocate some space
    counter = 1;
    %precisionLocs = {};
    %intensityLocs = {};
    %% perform filter in XYZ with either polymask or abs values
    %polymask
    if options.XYZ.XY.UseManualPolymask
        %perform cellmask
        %make folder to save the mask
        mkdir Cellmask;
        cd Cellmask;
        for i = 1:size(LocalisationData,1)
            tmpdata = LocalisationData{i,1};
            maxX = max(tmpdata(:,2));
            maxY = max(tmpdata(:,3));
            maxXY = max([maxX, maxY]);
            try
                tmpmaskdat = CellMaskV2(tmpdata, [round(maxXY*1.1) round(maxXY*1.1)], pwd, i);
            catch
                tmpmaskdat = {[1,1,1,1,1,1,1,1], [1,1,1,1,1,1,1,1]};
            end
            for j = 1:size(tmpmaskdat,1)
                maskedLocs{counter,1} = tmpmaskdat{j,1};
                newFileName = join([LocalisationData{i,2}, "_", int2str(j)],"");
                maskedLocs{counter,2} = newFileName;
                counter = counter+1;
            end
        end
        %remove leftover empty cells
        %maskedLocs{counter+1:end,:} = [];
        emptIdx = (cellfun('isempty',maskedLocs));
        maskedLocs(emptIdx) = [];
        lenLocs = length(maskedLocs);
        maskedLocs = reshape(maskedLocs,lenLocs/2,[]);

        cd ..;
        %if something went wrong, remove the file
        try
            for i = 1:size(maskedLocs,1)
                tmp = maskedLocs{i,1};
                if size(tmp,1)<5
                    maskedLocs(i,:) = [];
                end
            end
        catch
        end
        %absoltue values
    elseif options.XYZ.XY.X.Used || options.XYZ.XY.Y.Used
        %use the filter in XY
        for i = 1:size(LocalisationData,1)
            tmpdata = LocalisationData{i,1};
            if options.XYZ.XY.X.Used
                %filter by X value
                usrminx = options.XYZ.XY.X.Min;
                usrmaxx = options.XYZ.XY.X.Max;
                tmpdata = tmpdata(tmpdata(:,3)>= usrminx & tmpdata(:,3)<= usrmaxx,:);
            end
            if options.XYZ.XY.Y.Used
                %filter by Y value
                usrminy = options.XYZ.XY.Y.Min;
                usrmaxy = options.XYZ.XY.Y.Max;
                tmpdata = tmpdata(tmpdata(:,4)>= usrminy & tmpdata(:,4)<= usrmaxy,:);
            end
            maskedLocs{i,1} = tmpdata;
            maskedLocs{i,2} = LocalisationData{i,2};
        end
        %remove placeholders
        maskedLocs = maskedLocs(~cellfun(@isempty, maskedLocs));
        lenLocs = length(maskedLocs);
        maskedLocs = reshape(maskedLocs,lenLocs/2,[]);
        try
            for i = 1:size(maskedLocs,1)
                tmp = maskedLocs{i,1};
                if size(tmp,1)<5
                    maskedLocs(i,:) = [];
                end
            end
        catch
        end
    else
        %neither XY or polymask
        maskedLocs = LocalisationData;
    end
    %then do the Z
    if options.XYZ.Z.Used
        for i = 1:size(maskedLocs,1)
            tmpdata = maskedLocs{i,1};
            usrminz = options.XYZ.Z.Min;
            usrmaxz = options.XYZ.Z.Max;
            tmpdata = tmpdata(tmpdata(:,4)>= usrminz & tmpdata(:,4)<= usrmaxz,:);
            maskedLocs{i,1} = tmpdata;
        end
        try
            for i = 1:size(maskedLocs,1)
                tmp = maskedLocs{i,1};
                if size(tmp,1)<5
                    maskedLocs(i,:) = [];
                end
            end
        catch
        end
    end
    %% filter by precisions
    %XY precision
    precisionLocs = cell(size(maskedLocs));
    if options.precision.XY.Used
        for i = 1:size(maskedLocs,1)
            try
                tmpdata = maskedLocs{i,1};
                usrsigma = options.precision.XY.sigma;
                tmpdata = tmpdata(tmpdata(:,6)<= usrsigma & tmpdata(:,6)<= usrsigma,:);
                precisionLocs{i,1} = tmpdata;
                precisionLocs{i,2} = maskedLocs{i,2};
            catch
            end
        end
    else
        precisionLocs = maskedLocs;
    end
    %if something went wrong, remove the file
    try    
        for i = 1:size(precisionLocs,1)
                tmp = precisionLocs{i,1};
                if size(tmp,1)<5
                    precisionLocs(i,:) = [];
                end
        end
    catch
    end
    %Z precision
    if options.precision.Z.Used
        for i = 1:size(precisionLocs,1)
            try
            tmpdata = precisionLocs{i,1};
            usrsigma = options.precision.Z.sigma;
            tmpdata = tmpdata(tmpdata(:,7)<= usrsigma,:);
            precisionLocs{i,1} = tmpdata;
            precisionLocs{i,2} = precisionLocs{i,2};
            catch
            end
        end
    end
    %if something went wrong, remove the file
        try    
            for i = 1:size(precisionLocs,1)
                tmp = precisionLocs{i,1};
                if size(tmp,1)<5
                    precisionLocs(i,:) = [];
                end
            end
        catch
        end
    if ~options.precision.XY.Used && ~options.precision.Z.Used
        precisionLocs = maskedLocs;
    end
    %% filter by intensity
    intensityLocs = cell(size(precisionLocs));
    if options.intensity.Used
        for i = 1:size(precisionLocs,1)
            try
            tmpdata = precisionLocs{i,1};
            usrminInt = options.intensity.Min;
            usrmaxInt = options.intensity.Max; 
            tmpdata = tmpdata(tmpdata(:,8)>= usrminInt & tmpdata(:,8)<= usrmaxInt,:);
            intensityLocs{i,1} = tmpdata;
            intensityLocs{i,2} = precisionLocs{i,2};
            catch
            end
        end
        try
            for i = 1:size(intensityLocs,1)
                tmp = intensityLocs{i,1};
                if size(tmp,1)<5
                    intensityLocs(i,:) = [];
                end
            end
        catch
        end
    else
        intensityLocs = precisionLocs;
    end
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
    clear tmpdata tmpmaskdat;
    %save the localisations in mat and for swift
    datatosave = LocalisationData;
    save("SMAP_Localisations.mat","datatosave");

    datatosave = maskedLocs;
    save("maskedLocalisations.mat","datatosave");
    if saveCSV
        mkdir maskedLocsCSV
        tmpPath = pwd + "/maskedLocsCSV";
        for f = 1:size(datatosave,1)
            tmpFileData = datatosave{f,1};
            [~, tmpFileName, ~] = fileparts(datatosave{f,2});
            writematrix(tmpFileData, fullfile(tmpPath, tmpFileName, ".csv"));
        end
    end
    datatosave = precisionLocs;
    save("precisionFilteredLocalisations.mat","datatosave");
     if saveCSV
        mkdir precisionFilteredCSV
        tmpPath = pwd + "/precisionFilteredCSV";
        for f = 1:size(datatosave,1)
            tmpFileData = datatosave{f,1};
            [~, tmpFileName, ~] = fileparts(datatosave{f,2});
            writematrix(tmpFileData, fullfile(tmpPath, tmpFileName, ".csv"));
        end
    end
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