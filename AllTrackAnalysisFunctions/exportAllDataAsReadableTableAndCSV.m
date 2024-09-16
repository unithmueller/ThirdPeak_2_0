function exportAllDataAsReadableTableAndCSV(Trackdata, UnitCalculatedDataStruct, PixelCalculatedDataStruct)
%function to save the calculated values as readable csv files and tables to
%be used outside the software

%% Trackdata
curPath = pwd();
mkdir inputTrackData
addPath = "inputTrackData";

for i = 1:size(Trackdata,1)
    tmpName = Trackdata{i,2};
    tmpData = Trackdata{i,1};
    datSize = size(tmpData,2);
    if datSize <= 5 %minimal amount of data
        varNames = {"ID", "T", "X", "Y", "Z"};
    elseif datSize > 5 & datSize <= 23 & mean(tmpData(:,6:22),"all") ~= 0 %from swift/external tracking
        varNames = {"ID", "T", "X", "Y", "Z", "JumpDIst", "SegmentD", "SegmentDErr", "SegementDR", "SegmentMeanJumpDist", "SegmentMeanJumpDistErr", "SegmentMSD", "SegmentMSDErr", "Motiontype", "SegmentStart", "SegmentLifetime", "NumberOfSegments", "Intesity", "SegmentMeanJumpDIst", "Timelag", "TrackID", "Tracksegments"};
        varNames = varNames(1:datSize);
    elseif datSize > 23 & mean(tmpData(:,6:22),"all") == 0 %only DeepSPT
        varNames = {"ID", "T", "X", "Y", "Z", 'Alpha', 'D', 'extra', 'pval', 'Efficiency', 'logEfficiency', 'FractalDim','Gaussianity', 'Kurtosis', 'MSDratio', 'Trappedness', 't0', 't1', 't2', 't3', 'lifetime',  'length', 'avgSL', 'avgMSD', 'AvgDP', 'corrDP','signDP', 'sumSL', 'minSL', 'maxSL', 'BroadnessSL', 'Speed', 'CoV', 'FractionSlow',  'FractionFast', 'Volume', 'perc_ND', 'perc_DM', 'perc_CD', 'perc_SD', 'num_changepoints', 'inst_msd_D', 'meanSequence', 'medianSequence', 'maxSequence',  'minSequence', 'stdSequence', 'simSeq'};
    elseif datSize > 66 %external + DeepSPT
        varNames = {"ID", "T", "X", "Y", "Z", "JumpDIst", "SegmentD", "SegmentDErr", "SegementDR", "SegmentMeanJumpDist", "SegmentMeanJumpDistErr", "SegmentMSD", "SegmentMSDErr", "Motiontype", "SegmentStart", "SegmentLifetime", "NumberOfSegments", "Intesity", "SegmentMeanJumpDIst", "Timelag", "TrackID", "Tracksegments" 'Alpha', 'D', 'extra', 'pval', 'Efficiency', 'logEfficiency', 'FractalDim','Gaussianity', 'Kurtosis', 'MSDratio', 'Trappedness', 't0', 't1', 't2', 't3', 'lifetime',  'length', 'avgSL', 'avgMSD', 'AvgDP', 'corrDP','signDP', 'sumSL', 'minSL', 'maxSL', 'BroadnessSL', 'Speed', 'CoV', 'FractionSlow',  'FractionFast', 'Volume', 'perc_ND', 'perc_DM', 'perc_CD', 'perc_SD', 'num_changepoints', 'inst_msd_D', 'meanSequence', 'medianSequence', 'maxSequence',  'minSequence', 'stdSequence', 'simSeq', "segID"};

    end
    varNames = cellfun(@convertStringsToChars, varNames, 'UniformOutput', false);
    tmpData = array2table(tmpData,"VariableNames",varNames);
    [~,filename,~] = fileparts(tmpName);
    filename = filename+".csv";
    compFilename = fullfile(curPath, addPath, filename);
    writetable(tmpData, compFilename);
end
clear tmpData Trackdata
%% data struct
mkdir calculatedTrackStat
addPath = "calculatedTrackStat";

totalStrucDataSize = sum(cellfun(@(x) size(x,1), UnitCalculatedDataStruct.JumpDist.X(:,2)));

varNames = {"ID", "T", "JDX", "JDY", "JDZ", "JDXY", "JDXYZ", "CumJDX", "CumJDY", "CumJDZ", "CumJDXY", "CumJDXYZ", "MeanJDX", "MeanJDY", "MeanJDZ", "MeanJDXY", "MeanJDXYZ", "CumMeanJDX", "CumMeanJDY", "CumMeanJDZ", "CumMeanJDXY", "CumMeanJDXYZ", "JumpAnglesXY", "JumpAnglesXYZ", "swiftParMSD", "swiftParMSDerr", "swiftParD", "swiftParDerr", "swiftParType", "swiftSwitchFreq", "NumbrSteps", "NetLenXY", "NetLenXYZ", "ConfRatioXY", "ConfRatioXYZ", "AbsLenXY", "AbsLenXYZ", "MSDXYAlpha", "MSDXYA", "MSDXYd", "MSDXYLogR", "MSDXYLinR", "MSDXYZAlpha", "MSDXYZA", "MSDXYZd", "MSDXYZLogR", "MSDXYZLinR"};
varNames = cellfun(@convertStringsToChars, varNames, 'UniformOutput', false);

saveUnitArray = zeros(totalStrucDataSize, 47);
savePixelArray = zeros(totalStrucDataSize, 47);

%jumpDistance
JDXSZ = cellfun(@(x) size(x,1), UnitCalculatedDataStruct.JumpDist.X(:,2));
JDdat = UnitCalculatedDataStruct.JumpDist.X(:,2);
IDs = UnitCalculatedDataStruct.JumpDist.X(:,1);
TPdat = vertcat(JDdat{:});

idArray = zeros(size(TPdat,1),1);
%X
counter = 1;
for j = 1:size(IDs,1)
    tmpID = IDs{j,1};
    tmpSZ = JDXSZ(1);
    tmpIDArray = ones(tmpSZ,1)*tmpID;
    idArray(counter:counter+tmpSZ-1) = tmpIDArray;
    counter = counter+tmpSZ;
end
%id
saveUnitArray(:,1) = idArray;
savePixelArray(:,1) = idArray;
%TX
JDdat = PixelCalculatedDataStruct.JumpDist.X(:,2);
TPdat = vertcat(JDdat{:});
savePixelArray(:,2:3) = TPdat;

JDdat = UnitCalculatedDataStruct.JumpDist.X(:,2);
TPdat = vertcat(JDdat{:});
saveUnitArray(:,2:3) = TPdat;
%Y
JDdat = PixelCalculatedDataStruct.JumpDist.Y(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
savePixelArray(:,4) = JDdat;

JDdat = UnitCalculatedDataStruct.JumpDist.Y(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
saveUnitArray(:,4) = JDdat;
%Z
JDdat = PixelCalculatedDataStruct.JumpDist.Z(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
savePixelArray(:,5) = JDdat;

JDdat = UnitCalculatedDataStruct.JumpDist.Z(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
saveUnitArray(:,5) = JDdat;
%XY
JDdat = PixelCalculatedDataStruct.JumpDist.XY(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
savePixelArray(:,6) = JDdat;

JDdat = UnitCalculatedDataStruct.JumpDist.XY(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
saveUnitArray(:,6) = JDdat;
%XYZ
JDdat = PixelCalculatedDataStruct.JumpDist.XYZ(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
savePixelArray(:,7) = JDdat;

JDdat = UnitCalculatedDataStruct.JumpDist.XYZ(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,2);
saveUnitArray(:,7) = JDdat;
%CumJD
%X
JDdat = PixelCalculatedDataStruct.CumJumpDist.X(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
savePixelArray(:,8) = TPdat;

JDdat = UnitCalculatedDataStruct.CumJumpDist.X(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,8) = TPdat;
%Y
JDdat = PixelCalculatedDataStruct.CumJumpDist.Y(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
savePixelArray(:,9) = TPdat;

JDdat = UnitCalculatedDataStruct.CumJumpDist.Y(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,9) = TPdat;
%Z
JDdat = PixelCalculatedDataStruct.CumJumpDist.Z(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
savePixelArray(:,10) = TPdat;

JDdat = UnitCalculatedDataStruct.CumJumpDist.Z(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,10) = TPdat;
%XY
JDdat = PixelCalculatedDataStruct.CumJumpDist.XY(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
savePixelArray(:,11) = TPdat;

JDdat = UnitCalculatedDataStruct.CumJumpDist.XY(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,11) = TPdat;
%XYZ
JDdat = PixelCalculatedDataStruct.CumJumpDist.XYZ(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
savePixelArray(:,12) = TPdat;

JDdat = UnitCalculatedDataStruct.CumJumpDist.XYZ(:,2);
TPdat = vertcat(JDdat{:});
JDdat = TPdat(:,1);
saveUnitArray(:,12) = TPdat;

%meanJumpDist
%X Y Z XY XYZ
meanJDPixel = zeros(totalStrucDataSize,6);
meanJDPixel(:,1) = saveUnitArray(:,1);
meanJDUnit = zeros(totalStrucDataSize,6);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = PixelCalculatedDataStruct.MeanJumpDist.X{1,1};
pixlXDat = PixelCalculatedDataStruct.MeanJumpDist.X{1,2};
pixlYDat = PixelCalculatedDataStruct.MeanJumpDist.Y{1,2};
pixlZDat = PixelCalculatedDataStruct.MeanJumpDist.Z{1,2};
pixlXYDat = PixelCalculatedDataStruct.MeanJumpDist.XY{1,2};
pixlXYZDat = PixelCalculatedDataStruct.MeanJumpDist.XYZ{1,2};
UnitXDat = UnitCalculatedDataStruct.MeanJumpDist.X{1,2};
UnitYDat = UnitCalculatedDataStruct.MeanJumpDist.Y{1,2};
UnitZDat = UnitCalculatedDataStruct.MeanJumpDist.Z{1,2};
UnitXYDat = UnitCalculatedDataStruct.MeanJumpDist.XY{1,2};
UnitXYZDat = UnitCalculatedDataStruct.MeanJumpDist.XYZ{1,2};

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};
    meanJDPixel(meanJDPixel(:,1)==idx,2) = pixlXDat(j);
    meanJDPixel(meanJDPixel(:,1)==idx,3) = pixlYDat(j);
    meanJDPixel(meanJDPixel(:,1)==idx,4) = pixlZDat(j);
    meanJDPixel(meanJDPixel(:,1)==idx,5) = pixlXYDat(j);
    meanJDPixel(meanJDPixel(:,1)==idx,6) = pixlXYZDat(j);

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat(j);
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat(j);
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat(j);
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat(j);
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat(j);
end

saveUnitArray(:,13:17) = meanJDUnit(:,2:6);
savePixelArray(:,13:17) = meanJDPixel(:,2:6);

%cumMeanJD
%X Y Z XY XYZ
meanJDPixel = zeros(totalStrucDataSize,6);
meanJDPixel(:,1) = saveUnitArray(:,1);
meanJDUnit = zeros(totalStrucDataSize,6);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = PixelCalculatedDataStruct.MeanJumpDist.X{1,1};
pixlXDat = PixelCalculatedDataStruct.CumMeanJumpDist.X(:,2);
pixlYDat = PixelCalculatedDataStruct.CumMeanJumpDist.Y(:,2);
pixlZDat = PixelCalculatedDataStruct.CumMeanJumpDist.Z(:,2);
pixlXYDat = PixelCalculatedDataStruct.CumMeanJumpDist.XY(:,2);
pixlXYZDat = PixelCalculatedDataStruct.CumMeanJumpDist.XYZ(:,2);
UnitXDat = UnitCalculatedDataStruct.CumMeanJumpDist.X(:,2);
UnitYDat = UnitCalculatedDataStruct.CumMeanJumpDist.Y(:,2);
UnitZDat = UnitCalculatedDataStruct.CumMeanJumpDist.Z(:,2);
UnitXYDat = UnitCalculatedDataStruct.CumMeanJumpDist.XY(:,2);
UnitXYZDat = UnitCalculatedDataStruct.CumMeanJumpDist.XYZ(:,2);

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};
    meanJDPixel(meanJDPixel(:,1)==idx,2) = pixlXDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,3) = pixlYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,4) = pixlZDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,5) = pixlXYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,6) = pixlXYZDat{j};

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
end

saveUnitArray(:,18:22) = meanJDUnit(:,2:6);
savePixelArray(:,18:22) = meanJDPixel(:,2:6);

%jump angles
%XY
JDdat = PixelCalculatedDataStruct.JumpAngles.XY(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
savePixelArray(:,23) = TPdat;

JDdat = UnitCalculatedDataStruct.JumpAngles.XY(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
saveUnitArray(:,23) = TPdat;

%XYZ
JDdat = PixelCalculatedDataStruct.JumpAngles.XYZ(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
savePixelArray(:,24) = TPdat;

JDdat = UnitCalculatedDataStruct.JumpAngles.XYZ(:,2);
for j = 1:size(JDdat,1)
    tmp = JDdat{j,1};
    tmp(end+1) = 0;
    JDdat{j,1} = tmp;
end
TPdat = vertcat(JDdat{:});
saveUnitArray(:,24) = TPdat;

%swiftParams
%MSD MErr D Derr Type Freq
meanJDPixel = zeros(totalStrucDataSize,7);
meanJDPixel(:,1) = saveUnitArray(:,1);
meanJDUnit = zeros(totalStrucDataSize,7);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = PixelCalculatedDataStruct.MeanJumpDist.X{1,1};
pixlXDat = PixelCalculatedDataStruct.SwiftParams.MSD(:,2);
pixlYDat = PixelCalculatedDataStruct.SwiftParams.MSDerr(:,2);
pixlZDat = PixelCalculatedDataStruct.SwiftParams.D(:,2);
pixlXYDat = PixelCalculatedDataStruct.SwiftParams.Derr(:,2);
pixlXYZDat = PixelCalculatedDataStruct.SwiftParams.type(:,2);
pixlswitchFrq = PixelCalculatedDataStruct.SwiftParams.switchFreq(:,2);

UnitXDat = UnitCalculatedDataStruct.SwiftParams.MSD(:,2);
UnitYDat = UnitCalculatedDataStruct.SwiftParams.MSDerr(:,2);
UnitZDat = UnitCalculatedDataStruct.SwiftParams.D(:,2);
UnitXYDat = UnitCalculatedDataStruct.SwiftParams.Derr(:,2);
UnitXYZDat = UnitCalculatedDataStruct.SwiftParams.type(:,2);
Unitswitchfreq = UnitCalculatedDataStruct.SwiftParams.switchFreq(:,2);

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};
    meanJDPixel(meanJDPixel(:,1)==idx,2) = pixlXDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,3) = pixlYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,4) = pixlZDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,5) = pixlXYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,6) = pixlXYZDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,7) = pixlswitchFrq{j};

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,7) = Unitswitchfreq{j};
end

saveUnitArray(:,25:30) = meanJDUnit(:,2:7);
savePixelArray(:,25:30) = meanJDPixel(:,2:7);

%stepNumber
meanJDPixel = zeros(totalStrucDataSize,8);
meanJDPixel(:,1) = saveUnitArray(:,1);
meanJDUnit = zeros(totalStrucDataSize,8);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = PixelCalculatedDataStruct.MeanJumpDist.X{1,1};
pixlXDat = PixelCalculatedDataStruct.TrackLength.Steps(:,2);
pixlYDat = PixelCalculatedDataStruct.TrackLength.NetLength.XY(:,2);
pixlZDat = PixelCalculatedDataStruct.TrackLength.NetLength.XYZ(:,2);
pixlXYDat = PixelCalculatedDataStruct.TrackLength.ConfRatio.XY(:,2);
pixlXYZDat = PixelCalculatedDataStruct.TrackLength.ConfRatio.XYZ(:,2);
pixlabsLengthxy = PixelCalculatedDataStruct.TrackLength.AbsLength.XY(:,2);
pixlabsLengthxyz = PixelCalculatedDataStruct.TrackLength.AbsLength.XYZ(:,2);

UnitXDat = UnitCalculatedDataStruct.TrackLength.Steps(:,2);
UnitYDat = UnitCalculatedDataStruct.TrackLength.NetLength.XY(:,2);
UnitZDat = UnitCalculatedDataStruct.TrackLength.NetLength.XYZ(:,2);
UnitXYDat = UnitCalculatedDataStruct.TrackLength.ConfRatio.XY(:,2);
UnitXYZDat = UnitCalculatedDataStruct.TrackLength.ConfRatio.XYZ(:,2);
unitabsLengthxy = UnitCalculatedDataStruct.TrackLength.AbsLength.XY(:,2);
unitabsLengthxyz = UnitCalculatedDataStruct.TrackLength.AbsLength.XYZ(:,2);

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};
    meanJDPixel(meanJDPixel(:,1)==idx,2) = pixlXDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,3) = pixlYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,4) = pixlZDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,5) = pixlXYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,6) = pixlXYZDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,7) = pixlabsLengthxy{j};
    meanJDPixel(meanJDPixel(:,1)==idx,8) = pixlabsLengthxyz{j};

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,7) = unitabsLengthxy{j};
    meanJDUnit(meanJDUnit(:,1)==idx,8) = unitabsLengthxyz{j};
end

saveUnitArray(:,31:37) = meanJDUnit(:,2:8);
savePixelArray(:,31:37) = meanJDPixel(:,2:8);

%internMSDFit
%XY
meanJDPixel = zeros(totalStrucDataSize,6);
meanJDPixel(:,1) = saveUnitArray(:,1);
meanJDUnit = zeros(totalStrucDataSize,6);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = PixelCalculatedDataStruct.MeanJumpDist.X{1,1};
pixlXDat = PixelCalculatedDataStruct.InternMSD.XY.Alpha(:,2);
pixlYDat = PixelCalculatedDataStruct.InternMSD.XY.a(:,2);
pixlZDat = PixelCalculatedDataStruct.InternMSD.XY.d(:,2);
pixlXYDat = PixelCalculatedDataStruct.InternMSD.XY.logR(:,2);
pixlXYZDat = PixelCalculatedDataStruct.InternMSD.XY.linR(:,2);

UnitXDat = UnitCalculatedDataStruct.InternMSD.XY.Alpha(:,2);
UnitYDat = UnitCalculatedDataStruct.InternMSD.XY.a(:,2);
UnitZDat = UnitCalculatedDataStruct.InternMSD.XY.d(:,2);
UnitXYDat = UnitCalculatedDataStruct.InternMSD.XY.logR(:,2);
UnitXYZDat = UnitCalculatedDataStruct.InternMSD.XY.linR(:,2);

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};
    meanJDPixel(meanJDPixel(:,1)==idx,2) = pixlXDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,3) = pixlYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,4) = pixlZDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,5) = pixlXYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,6) = pixlXYZDat{j};

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
end
saveUnitArray(:,38:42) = meanJDUnit(:,2:6);
savePixelArray(:,38:42) = meanJDPixel(:,2:6);

%XYZ
meanJDPixel = zeros(totalStrucDataSize,6);
meanJDPixel(:,1) = saveUnitArray(:,1);
meanJDUnit = zeros(totalStrucDataSize,6);
meanJDUnit(:,1) = saveUnitArray(:,1);
IDXs = PixelCalculatedDataStruct.MeanJumpDist.X{1,1};
pixlXDat = PixelCalculatedDataStruct.InternMSD.XYZ.Alpha(:,2);
pixlYDat = PixelCalculatedDataStruct.InternMSD.XYZ.a(:,2);
pixlZDat = PixelCalculatedDataStruct.InternMSD.XYZ.d(:,2);
pixlXYDat = PixelCalculatedDataStruct.InternMSD.XYZ.logR(:,2);
pixlXYZDat = PixelCalculatedDataStruct.InternMSD.XYZ.linR(:,2);

UnitXDat = UnitCalculatedDataStruct.InternMSD.XYZ.Alpha(:,2);
UnitYDat = UnitCalculatedDataStruct.InternMSD.XYZ.a(:,2);
UnitZDat = UnitCalculatedDataStruct.InternMSD.XYZ.d(:,2);
UnitXYDat = UnitCalculatedDataStruct.InternMSD.XYZ.logR(:,2);
UnitXYZDat = UnitCalculatedDataStruct.InternMSD.XYZ.linR(:,2);

for j = 1:size(IDXs,1)
    idx = IDXs{j,1};
    meanJDPixel(meanJDPixel(:,1)==idx,2) = pixlXDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,3) = pixlYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,4) = pixlZDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,5) = pixlXYDat{j};
    meanJDPixel(meanJDPixel(:,1)==idx,6) = pixlXYZDat{j};

    meanJDUnit(meanJDUnit(:,1)==idx,2) = UnitXDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,3) = UnitYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,4) = UnitZDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,5) = UnitXYDat{j};
    meanJDUnit(meanJDUnit(:,1)==idx,6) = UnitXYZDat{j};
end
saveUnitArray(:,43:47) = meanJDUnit(:,2:6);
savePixelArray(:,43:47) = meanJDPixel(:,2:6);

%% save the table
saveUnitArray = array2table(saveUnitArray,"VariableNames",varNames);
savePixelArray = array2table(savePixelArray,"VariableNames",varNames);

UnitName = "UnitSaveStructure.csv";
pixelName = "PixelSaveStructure.csv";

compFilename = fullfile(curPath, addPath, UnitName);
writetable(saveUnitArray, compFilename);
compFilename = fullfile(curPath, addPath, pixelName);
writetable(savePixelArray, compFilename);

end