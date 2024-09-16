function flatUnitData = generateFlatOwnCalculatedData(Trackdata, UnitCalculatedDataStruct)
%flatten the structered array values to be used for supermap generation in
%the main window
Trackdata = Trackdata{1,1};

totalStrucDataSize = size(Trackdata,1);

varNames = {"ID", "T", "X", "Y", "Z", "JDX", "JDY", "JDZ", "JDXY", "JDXYZ", "CumJDX", "CumJDY", "CumJDZ", "CumJDXY", "CumJDXYZ", "MeanJDX", "MeanJDY", "MeanJDZ", "MeanJDXY", "MeanJDXYZ", "CumMeanJDX", "CumMeanJDY", "CumMeanJDZ", "CumMeanJDXY", "CumMeanJDXYZ", "JumpAnglesXY", "JumpAnglesXYZ", "NumbrSteps", "NetLenXY", "NetLenXYZ", "ConfRatioXY", "ConfRatioXYZ", "AbsLenXY", "AbsLenXYZ", "MSDXYAlpha", "MSDXYA", "MSDXYd", "MSDXYLogR", "MSDXYLinR", "MSDXYZAlpha", "MSDXYZA", "MSDXYZd", "MSDXYZLogR", "MSDXYZLinR"};
varNames = cellfun(@convertStringsToChars, varNames, 'UniformOutput', false);

saveUnitArray = zeros(totalStrucDataSize, 45);
%id, t, coordinates
saveUnitArray(:,1:5) = Trackdata(:,1:5);
%jumpDistance
JDXSZ = cellfun(@(x) size(x,1), UnitCalculatedDataStruct.JumpDist.X(:,2));
JDdat = UnitCalculatedDataStruct.JumpDist.X(:,2);
IDs = UnitCalculatedDataStruct.JumpDist.X(:,1);
TPdat = vertcat(JDdat{:});

%% JD
%X
JDdat = UnitCalculatedDataStruct.JumpDist.X(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1,:) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
saveUnitArray(:,6) = TPdat(:,2);
%Y
JDdat = UnitCalculatedDataStruct.JumpDist.Y(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1,:) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
saveUnitArray(:,7) = JDdat;
%Z
JDdat = UnitCalculatedDataStruct.JumpDist.Z(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1,:) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
saveUnitArray(:,8) = JDdat;
%XY
JDdat = UnitCalculatedDataStruct.JumpDist.XY(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1,:) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
saveUnitArray(:,9) = JDdat;
%XYZ
JDdat = UnitCalculatedDataStruct.JumpDist.XYZ(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1,:) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
saveUnitArray(:,10) = JDdat;
%% CumJD
%X
JDdat = UnitCalculatedDataStruct.CumJumpDist.X(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,11) = TPdat;
%Y
JDdat = UnitCalculatedDataStruct.CumJumpDist.Y(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,12) = TPdat;
%Z
JDdat = UnitCalculatedDataStruct.CumJumpDist.Z(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,13) = TPdat;
%XY
JDdat = UnitCalculatedDataStruct.CumJumpDist.XY(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,14) = TPdat;
%XYZ
JDdat = UnitCalculatedDataStruct.CumJumpDist.XYZ(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,15) = TPdat;

%% meanJumpDist
%X Y Z XY XYZ

meanJDUnit = zeros(totalStrucDataSize,6);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = UnitCalculatedDataStruct.MeanJumpDist.X{1,1};

UnitXDat = UnitCalculatedDataStruct.MeanJumpDist.X{1,2};
UnitYDat = UnitCalculatedDataStruct.MeanJumpDist.Y{1,2};
UnitZDat = UnitCalculatedDataStruct.MeanJumpDist.Z{1,2};
UnitXYDat = UnitCalculatedDataStruct.MeanJumpDist.XY{1,2};
UnitXYZDat = UnitCalculatedDataStruct.MeanJumpDist.XYZ{1,2};

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat(j);
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat(j);
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat(j);
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat(j);
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat(j);
end

saveUnitArray(:,16:20) = meanJDUnit(:,2:6);

%% cumMeanJD
%X Y Z XY XYZ

meanJDUnit = zeros(totalStrucDataSize,6);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = UnitCalculatedDataStruct.MeanJumpDist.X{1,1};

UnitXDat = UnitCalculatedDataStruct.CumMeanJumpDist.X(:,2);
UnitYDat = UnitCalculatedDataStruct.CumMeanJumpDist.Y(:,2);
UnitZDat = UnitCalculatedDataStruct.CumMeanJumpDist.Z(:,2);
UnitXYDat = UnitCalculatedDataStruct.CumMeanJumpDist.XY(:,2);
UnitXYZDat = UnitCalculatedDataStruct.CumMeanJumpDist.XYZ(:,2);

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
end

saveUnitArray(:,21:25) = meanJDUnit(:,2:6);

%% jump angles
%XY

JDdat = UnitCalculatedDataStruct.JumpAngles.XY(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1:end+2) = [0;0];
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
saveUnitArray(:,26) = TPdat;

%XYZ

JDdat = UnitCalculatedDataStruct.JumpAngles.XYZ(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1:end+2) = [0;0];
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
saveUnitArray(:,27) = TPdat;

%% stepNumber

meanJDUnit = zeros(totalStrucDataSize,8);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = UnitCalculatedDataStruct.MeanJumpDist.X{1,1};

UnitXDat = UnitCalculatedDataStruct.TrackLength.Steps(:,2);
UnitYDat = UnitCalculatedDataStruct.TrackLength.NetLength.XY(:,2);
UnitZDat = UnitCalculatedDataStruct.TrackLength.NetLength.XYZ(:,2);
UnitXYDat = UnitCalculatedDataStruct.TrackLength.ConfRatio.XY(:,2);
UnitXYZDat = UnitCalculatedDataStruct.TrackLength.ConfRatio.XYZ(:,2);
unitabsLengthxy = UnitCalculatedDataStruct.TrackLength.AbsLength.XY(:,2);
unitabsLengthxyz = UnitCalculatedDataStruct.TrackLength.AbsLength.XYZ(:,2);

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,7) = unitabsLengthxy{j};
    meanJDUnit(meanJDUnit(:,1)==idx,8) = unitabsLengthxyz{j};
end

saveUnitArray(:,28:34) = meanJDUnit(:,2:8);


%% internMSDFit
%XY
try
    meanJDUnit = zeros(totalStrucDataSize,6);
    meanJDUnit(:,1) = saveUnitArray(:,1);
    IDXs = UnitCalculatedDataStruct.MeanJumpDist.X{1,1};
    
    
    UnitXDat = UnitCalculatedDataStruct.InternMSD.XY.Alpha(:,2);
    UnitYDat = UnitCalculatedDataStruct.InternMSD.XY.a(:,2);
    UnitZDat = UnitCalculatedDataStruct.InternMSD.XY.d(:,2);
    UnitXYDat = UnitCalculatedDataStruct.InternMSD.XY.logR(:,2);
    UnitXYZDat = UnitCalculatedDataStruct.InternMSD.XY.linR(:,2);
    
    for j = 1:size(IDXs,1)
        idx = IDXs{j,1};
    
    
        meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
        meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
        meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
        meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
        meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
    end
    saveUnitArray(:,35:40) = meanJDUnit(:,2:6);
    
    
    %XYZ
    
    meanJDUnit = zeros(totalStrucDataSize,6);
    meanJDUnit(:,1) = saveUnitArray(:,1);
    IDXs = UnitCalculatedDataStruct.MeanJumpDist.X{1,1};
    
    
    UnitXDat = UnitCalculatedDataStruct.InternMSD.XYZ.Alpha(:,2);
    UnitYDat = UnitCalculatedDataStruct.InternMSD.XYZ.a(:,2);
    UnitZDat = UnitCalculatedDataStruct.InternMSD.XYZ.d(:,2);
    UnitXYDat = UnitCalculatedDataStruct.InternMSD.XYZ.logR(:,2);
    UnitXYZDat = UnitCalculatedDataStruct.InternMSD.XYZ.linR(:,2);
    
    for j = 1:size(IDXs,1)
        idx = IDXs{j,1};
    
    
        meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
        meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
        meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
        meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
        meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
    end
    saveUnitArray(:,41:45) = meanJDUnit(:,2:6);
catch
end

flatUnitData = saveUnitArray;
end