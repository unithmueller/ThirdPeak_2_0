function [minv, maxv, gaussDat, kernelDat,newBinNmbrs] = plotCumulativeJumpDistance(Axes, binNumbers, SaveStructure, dimension, isPixel, lengthUnit, filterIDs, performFit, autoBins)
%Function to plot the respective cumulative jump distances as a histogram in the track
%analysis window.
%Input: Axes - axes object to plot into
       %SaveStructure - structured array that contains the jump distances
       %dimension - defines the dimensionality to be plotted (X,Y,Z,XY,XYZ)
       %filterIDs - if we use a filter we can provide the ids to search for
       %here
       %performFit - will try to perform several fits to find a matching
       %distribution
       
       %% Get the data of choice
       data = SaveStructure.CumJumpDist.(dimension);
       
       %% Apply the filter if necessary
       if size(filterIDs,1)>0
           ids = data(:,1);
           ids = cell2mat(ids);
           idx = ismember(ids, filterIDs);
           filteredData = cell(size(idx,1),2);
           for i = 1:size(idx)
               if idx(i)
                  filteredData{i,1} = data{i,1};
                  filteredData{i,2} = data{i,2};
               end
           end
           filteredData = filteredData((~cellfun(@(x) isempty(x), filteredData(:,1))),:);
           data = filteredData;
       end
       %% Unpack the cell array
       data = cell2mat(data(:,2));
       
       %% Plot the data
       minv = min(data);
       maxv = max(data);
       if autoBins
           binNum = round(sqrt(size(data,1)));
           edges = linspace(minv, maxv, binNum);
           newBinNmbrs = binNum;
       else
           edges = linspace(minv, maxv, binNumbers);
           newBinNmbrs = binNumbers;
       end
       his = histogram(Axes, data, edges);
       hixMaxValue = max(his.Values);
       xlim(Axes, [minv maxv]);
       title(Axes, join(["Cumulative Jump Distance Distribution for " dimension],""));
       if isPixel
           xlabel(Axes, "Cumulative Jump Distance in [px]");
       else
           xlabel(Axes, sprintf("Cumulative Jump Distance in [%s]", lengthUnit));
       end
       ylabel(Axes, "Counts");
       
       %% Decide if we fit or not
       if performFit
           %% perform fit
           pdGauss = fitdist(data, "Normal");
           pdKernel = fitdist(data, "Kernel", "Width", []);

           %% generate matching data
           xFitData = minv:1:maxv;
           yGauss = pdf(pdGauss, xFitData);
           yKernel = pdf(pdKernel, xFitData);
           maxyGauss = max(yGauss);
           maxyKernel = max(yKernel);
           %scaling factor
           scalingFactorGauss = hixMaxValue/maxyGauss;
           scalingFactorKernel = hixMaxValue/maxyKernel;
           
           %scale the data
           yGauss = yGauss*scalingFactorGauss;
           yKernel = yKernel*scalingFactorKernel;

           %% plot
           axes(Axes);
           hold(Axes,"on")
           gp = plot(Axes, xFitData, yGauss, "--r");
           kp = plot(Axes, xFitData, yKernel, "k");
           legend(Axes, "Histogram", "GaussFit", "KernelFit");
           hold(Axes,"off")
           
           %% get the data from fit
           gaussNLL = negloglik(pdGauss);
           kernelNLL = negloglik(pdKernel);
           gaussDat = [median(pdGauss), mean(pdGauss), std(pdGauss), (std(pdGauss)/sqrt(size(data,1))), var(pdGauss), gaussNLL];
           kernelDat = [median(pdKernel), mean(pdKernel), std(pdKernel), (std(pdKernel)/sqrt(size(data,1))), var(pdKernel), kernelNLL];
       else
           gaussDat = [];
           kernelDat = [];
       end
end