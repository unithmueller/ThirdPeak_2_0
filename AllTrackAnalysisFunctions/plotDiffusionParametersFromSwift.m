function [minv, maxv, gaussDat, kernelDat, newBinNmbrs] = plotDiffusionParametersFromSwift(Axes, binNumbers, SaveStructure, property, lengthUnit, filterIDs, performFit, timeunit, autoBins)
%Function to plot the respective jump distances as a histogram in the track
%analysis window.
%Input: Axes - axes object to plot into
       %SaveStructure - structured array that contains the jump distances
       %property - defines if either D or MSD
       %filterIDs - if we use a filter we can provide the ids to search for
       %here
       %performFit - will try to perform several fits to find a matching
       %distribution
       
       %% Get the data of choice
       data = SaveStructure.SwiftParams.(property);
       
       %% Apply the filter if necessary
       if size(filterIDs,1)>0
           ids = data(:,1);
           ids = cell2mat(ids);
           mask = ismember(ids, filterIDs);
           filteredData = data(mask,:);
           data = filteredData;
       end

       %% Unpack the cell array
       data = cell2mat(data(:,2));
       
       %% Plot the data
       if mean(data,2) == 0
           return
       end
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
       title(Axes, join(["Distribution of " property " from Swift"],""));
       if property == "MSD"
           xlabel(Axes, "MSD in [µm²/s²]");
       else
           xlabel(Axes, sprintf("D in [µm²/%s]", timeunit));
       end
 
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
           gaussDat = [median(pdGauss), mean(pdGauss), std(pdGauss), (std(pdGauss)/sqrt(size(data,1))), var(pdGauss), negloglik(pdGauss)];
           kernelDat = [median(pdKernel), mean(pdKernel), std(pdKernel), (std(pdKernel)/sqrt(size(data,1))), var(pdKernel), negloglik(pdKernel)];
       else
           gaussDat = [];
           kernelDat = [];
       end
end