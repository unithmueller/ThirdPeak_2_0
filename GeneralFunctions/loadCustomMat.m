function tracks = loadCustomMat(file, ImportSettingsStruct)
%Function to allow for custom .mat localisation files. Data columns need to
%be given by the user in the import settings window.

    file = load(file);
    name = fieldnames(file);
    name = name{1,1};
    file = file.(name);
    type = ImportSettingsStruct.islocs;
    
    if type
    %newfile(:,1) = [];
    newfile(:,2) = file(:,ImportSettingsStruct.framenr); %t
    newfile(:,3) = file(:,ImportSettingsStruct.xpos); %x
    newfile(:,4) = file(:,ImportSettingsStruct.ypos); %y
    try
        newfile(:,5) = file(:,ImportSettingsStruct.zpos); %z
    catch
    end
    if size(newfile,2) < 5
        newfile(:,5) = 0; %z
    end
    try
        newfile(:,6) = file(:,ImportSettingsStruct.xyerr); %xyerr
    catch
       newfile(:,6) = 0;
    end
    try
        newfile(:,7) = file(:,ImportSettingsStruct.zerr); %zerr
    catch
        newfile(:,7) = 0;
    end
    try
        newfile(:,8) = file(:,ImportSettingsStruct.int); %photons
    catch
        newfile(:,8) = 0;
    end
    try
        newfile(:,9) = file(:,ImportSettingsStruct.interr); %photon error
    catch
        newfile(:,9) = 0;
    end
    newfile(:,10:22) = 0;
    tracks = newfile;
        return
    else
    newfile(:,1) = file(:,ImportSettingsStruct.trackid); %id
    newfile(:,2) = file(:,ImportSettingsStruct.framenr); %t
    newfile(:,3) = file(:,ImportSettingsStruct.xpos); %x
    newfile(:,4) = file(:,ImportSettingsStruct.ypos); %y
    try
        newfile(:,5) = file(:,ImportSettingsStruct.zpos); %z
    catch
    end
    if size(newfile,2) < 5
        newfile(:,5) = 0; %z
    end
    %optional data, try to grab it. definetly not the most beautiful way of
    %implementing this....
    try
        newfile(:,6) = file(:,6); %jumpdist
    catch
        newfile(:,6) = 0;
    end
    try
        newfile(:,7) = file(:,app.ImportSettingsStruct.d); %segmentD
    catch
        newfile(:,7) = 0;
    end
    try
        newfile(:,8) = file(:,app.ImportSettingsStruct.derr); %segmentDerr
    catch
        newfile(:,8) = 0;
    end
    try
        newfile(:,9) = file(:,14); %segmentDR
    catch
        newfile(:,9) = 0;
    end
    try
        newfile(:,10) = file(:,15); %segmentmeanjumpdist
    catch
        newfile(:,10) = 0;
    end
    try
        newfile(:,11) = file(:,16); %segmeanjumpdisterr
    catch
        newfile(:,11) = 0;
    end
    try
        newfile(:,12) = file(:,app.ImportSettingsStruct.msd); %segment-MSD
    catch
        newfile(:,12) = 0;
    end
    try
        newfile(:,13) = file(:,app.ImportSettingsStruct.msderr); %segment_MSDERR
    catch
        newfile(:,13) = 0;
    end
    try
        difftypedata = file(:,ImportSettingsStruct.difftype);
        diffstates = file(:,ImportSettingsStruct.diffstates);
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
        newfile(:,15) = file(:,22); %segmentSTart
    catch
        newfile(:,15) = 0;
    end
    try
        newfile(:,16) = file(:,23); %segmentLifetime
    catch
        newfile(:,16) = 0;
    end
    try
        newfile(:,17) = file(:,26); %numberofsegmentintrack
    catch
        newfile(:,17) = 0;
    end
    newfile(:,18:22) = 0;
    %return the new tracks
    tracks = newfile;
    end
end