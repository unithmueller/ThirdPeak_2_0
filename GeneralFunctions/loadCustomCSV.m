function tracks = loadCustomCSV(file, ImportSettingsStruct)
%Function to load localisations from a custom CSV file. Import settings are
%given by the import settings struc by the user in the gui
    delimit = ImportSettingsStruct.separator;
    try
        delimit = convertStringsToChars(delimit);
    catch
    end
    skip = ImportSettingsStruct.headlerlines;
    readOpts = detectImportOptions(file);
    readOpts.Delimiter = {delimit};
    readOpts.DataLines = [skip+1, Inf];

    file = readtable(file, readOpts);
    %file = table2array(file);

    %file = dlmread(file, del, skip);
    
    type = ImportSettingsStruct.islocs;

    if type
    newfile(:,1) = zeros(size(file,1),1); %placeholder
    newfile(:,2) = table2array(file(:,ImportSettingsStruct.framenr)); %t
    newfile(:,3) = table2array(file(:,ImportSettingsStruct.xpos)); %x
    newfile(:,4) = table2array(file(:,ImportSettingsStruct.ypos)); %y
    try
        newfile(:,5) = table2array(file(:,ImportSettingsStruct.zpos)); %z
    catch
    end
    if size(newfile,2) < 5
        newfile(:,5) = 0;
    end
    try
        newfile(:,6) = table2array(file(:,ImportSettingsStruct.xyerr)); %xyerr
    catch
        newfile(:,6) = 0;
    end
    try
        newfile(:,7) = table2array(file(:,ImportSettingsStruct.zerr)); %zerr
    catch
        newfile(:,7) = 0;
    end
    try
        newfile(:,8) = table2array(file(:,ImportSettingsStruct.int)); %photons
    catch
        newfile(:,8) = 0;
    end
    try
        newfile(:,9) = table2array(file(:,ImportSettingsStruct.interr)); %photon error
    catch
        newfile(:,9) = 0;
    end
    try
        newfile(:,10) = table2array(file(:,ImportSettingsStruct.bg)); %background
    catch
        newfile(:,10) = 0;
    end
    newfile(:,11:22) = 0;
    tracks = newfile;
        return
    else
    newfile(:,1) = table2array(file(:,ImportSettingsStruct.trackid)); %id
    newfile(:,2) = table2array(file(:,ImportSettingsStruct.framenr)); %t
    newfile(:,3) = table2array(file(:,ImportSettingsStruct.xpos)); %x
    newfile(:,4) = table2array(file(:,ImportSettingsStruct.ypos)); %y
    try
        newfile(:,5) = table2array(file(:,ImportSettingsStruct.zpos)); %z
    catch
    end
    if size(newfile,2) < 5
        newfile(:,5) = 0;
    end
    try
        newfile(:,6) = table2array(file(:,ImportSettingsStruct.jumpdist)); %jumpdist
    catch
        newfile(:,6) = 0;
    end
    try
        newfile(:,7) = table2array(file(:,ImportSettingsStruct.d)); %segmentD
    catch
        newfile(:,7) = 0;
    end
    try
        newfile(:,8) = table2array(file(:,ImportSettingsStruct.derr)); %segmentDerr
    catch
        newfile(:,8) = 0;
    end
    try
        newfile(:,9) = table2array(file(:,14)); %segmentDR
    catch
        newfile(:,9) = 0;
    end
    try
        newfile(:,10) = table2array(file(:,ImportSettingsStruct.meanjumpdist)); %segmentmeanjumpdist
    catch
        newfile(:,10) = 0;
    end
    try
        newfile(:,11) = table2array(file(:,ImportSettingsStruct.meanjumpdisterr)); %segmeanjumpdisterr
    catch
        newfile(:,11) = 0;
    end
    try
        newfile(:,12) = table2array(file(:,ImportSettingsStruct.msd)); %segment-MSD
    catch
        newfile(:,12) = 0;
    end
    try
        newfile(:,13) = table2array(file(:,ImportSettingsStruct.msderr)); %segment_MSDERR
    catch
        newfile(:,13) = 0;
    end
    try
        difftypedata = table2array(file(:,ImportSettingsStruct.difftype));
        diffstates = table2array(file(:,ImportSettingsStruct.diffstates));
        numbrOfStates = size(diffstates,2);
        for i = 1:numbrOfStates
            searchterm = diffstates(i);
            difftypedata(difftypedata == searchterm) = i;
        end
        newfile(:,14) = difftypedata; %motiontype
    catch
        newfile(:,14) = 0;
    end
    try
        newfile(:,15) = table2array(file(:,22)); %segmentSTart
    catch
        newfile(:,15) = 0;
    end
    try
        newfile(:,16) = table2array(file(:,23)); %segmentLifetime
    catch
        newfile(:,16) = 0;
    end
    try
        newfile(:,17) = table2array(file(:,26)); %numberofsegmentintrack 
    catch
        newfile(:,17) = 0;
    end
    %return the new tracks
    newfile(:,18:22) = 0;
    if isa(newfile, "table")
        tracks = table2array(newfile);
    else
        tracks = newfile;
    end
end