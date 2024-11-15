function [minv, maxv, nBins gaussDat, kernelDat] = performJDDAnalysis(Axes, TrackData, dimension, isPixel, lengthUnit, pxsize, timeunit, timestep, filterIDs, isSlidingWindow, windowSize, binDetermineMode, selectedDiffusionModes, findBestModelIterations, threeDastwoD, doModelFitting, performFit)
%%Wrapping function to acess the JDD package and make it available for the
%%user via the grapical user interface.
a= 1;
if isempty(selectedDiffusionModes)
    return
end
    %% filter the data if necessary
    if size(filterIDs,1) > 0
        mask = ismember(TrackData(:,1), filterIDs);
        TrackData = TrackData(mask, :);
    end

    %% filter the data by the step size if track is shorter than the step
    trackIds = unique(TrackData(:,1));
    trackLengths = zeros(size(trackIds,1),2);
    for i = 1:size(trackIds,1)
        trackLengths(i,1) = trackIds(i);
        trackLengths(i,2) = size(TrackData(TrackData(:,1) == trackIds(i),1),1);
    end

    trackLengths = trackLengths(trackLengths(:,2) > windowSize+1,1);

    mask = ismember(TrackData(:,1), trackLengths);
    TrackData = TrackData(mask, :);

    %% Bring data in the right format, as well as deal with the slinding window and window size and calculate jump distances
    switch dimension
        case "X"
            dimension = 3;
        case "Y"
            dimension = 4;
        case "Z"
            dimension = 5;
        case "XY"
            dimension = [3,4];
        case "XYZ"
            dimension = [3,4,5];
    end

    [jumpDistances, x, y, z] = transformTrackDataIntoTimeWindows(TrackData, dimension, windowSize, isSlidingWindow, trackLengths,isPixel, timestep);
  
    %% Set the bin number
    %Number of Bins for fitting
    %Choose option here.
    N = size(jumpDistances,1);
    switch binDetermineMode
        case "SturgesRule"
            Nb=round(1+log2(N)); %Sturges Rule
        case "DoanesRule"
            sigma=sqrt(6*(N-2)/(N+1)/(N+3)); %for Doane's rule
            Nb=round(1+log2(N)+log2(1+abs(skewness(jumpDistances))/sigma)); %Doanes Rule
        case "RiceRule"
            Nb=round(2*(N^(1/3))); %Rice Rule
        case "SquareRootGuidance"
            Nb=round(sqrt(N)); %square root guidance
        case "ScottsNormalRule"
            Nb=round((max(jumpDistances)-min(jumpDistances))*N^(1/3)/(3.5*std(jumpDistances))); %Scott's Normal Reference Rule.
        case "FreedmanDiaconisRule"
            Nb=round((max(jumpDistances)-min(jumpDistances))*N^(1/3)/(2*iqr(jumpDistances))); %Freedman Diaconis Rule
    end

    %% Generate Jump distance histogram
    [dr, Ni, yi, ri]=BinningHist(jumpDistances, N, Nb,'no');
    
    %% Fit a suitable model
    points = windowSize;
    tau = (windowSize-1)*timestep;
    dt = timestep;

    [param, beta, betaBoot, modelProb, modelVal, modelMethod,strucNames, nameArray, indArray,diffModes] = fitAndDetermineJDDModel(selectedDiffusionModes, isSlidingWindow, findBestModelIterations, dimension, threeDastwoD, tau, dr, ri, yi, Ni, N, points, dt, jumpDistances,x,y,z, Nb, doModelFitting);

    %% Plot it

    [exportParamNameArray, exportParamArray, nBins, minv, maxv] = plotJDDHistogramAndEverythingElse(Axes, jumpDistances, Nb, dimension, param, beta, betaBoot, modelProb,strucNames, nameArray, indArray,diffModes,tau, dr, ri, isPixel, lengthUnit, pxsize, timeunit, threeDastwoD, N);
       
    %% Gather the data for the table if possible
    gaussDat = {exportParamArray, exportParamNameArray};
    try
        kernelDat = modelProb;
    catch
        kernelDat = [];
    end
end



