function processedTracks = processTracksByDeepSPT(trackfiles)
%% Need to send the data to python, however cells are not supported
trackfiles = testfiles;
numberOfFiles = length(trackfiles);
convertedFiles = py.list();
convertedPaths = py.list();

for i = 1:numberOfFiles
    tmp = py.numpy.array(trackfiles{i,1});
    tmpPath = trackfiles{i,2};
    convertedFiles.append(tmp);
    convertedPaths.append(tmpPath);
end

filePath = "C:\Users\auxin\Desktop\Pythontest";
dimension = 3;
timestep = 1/25;
windowSize = 5;
usedTmpSegNN = "";
usedHMM = "";
generateNew = 1;
%Packed the cell array components in a list with numpy arrays, need to
%unpack and restructure in python then
a = 1;

processedTracks = pyrunfile("ManageDeepSPTCallsfromMatlab.py", "results", matlabTracks=convertedFiles, fileNames=convertedPaths, filePath=filePath, dimension=dimension, timestep=timestep, windowSize=windowSize, usedTmpSegNN=usedTmpSegNN, usedHMM=usedHMM, generateNew=generateNew);

%% convert list of array to matrix
tmpResults = cell(processedTracks);
tmpResults = cellfun(@double, tmpResults, 'UniformOutput', false);
tmpResults = cell2mat(tmpResults);
tmpResults = reshape(tmpResults, [], 40);


%results = manageDeepSPTCallsfromMatlab(matlabTracks, fileNames, filePath, dimension, timestep, windowSize, usedTmpSegNN, usedHMM, generateNew)
%py.pickle.dump(convertedFiles, py.open("C:\\Users\\auxin\\Desktop\\Pythontest\\DeepSPT\\tracks.pkl", "wb"))

end