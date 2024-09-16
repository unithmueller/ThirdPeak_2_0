function [minv, maxv, gaussDat, kernelDat, newBinNmbrs] = plotDeepSPTDataAllTracks(Axes, binNumbers, Data, propertyID, propertynames, lengthUnit, timeUnit, filterIDs, performFit, filterByID, averageSwitch, autoBins, allPropIDs)
%function to manage the plotting of the deepspt data in the all track
%analysis window
%Input: Axes - axes object to plot into
       %bin numbers of the array to use
       %Data - data to work with
       %filterIDs - if we use a filter we can provide the ids to search for
       %here
       %length and time unit for the plotting
       %performFit - will try to perform several fits to find a matching
       %distribution
a = 1;

%% filter the data if necessary
    if size(filterIDs,1) > 0
        mask = ismember(Data(:,filterByID), filterIDs);
        Data = Data(mask, :);
    end
    
    %% grab the relevant property
    data = Data(:,[1,(propertyID)]);
    if averageSwitch
        newdata = zeros(size(unique(data(:,1)),1),2);
        newdata(:,1) = unique(data(:,1));
        for i = 1:size(unique(data(:,1)),1)
            tmpdat = data(data(:,1) == newdata(i,1), 2);
            tmpdat = mean(tmpdat);
            newdata(i,2) = tmpdat;
        end
        data = newdata;
    end
    %set the labels depending on the selected property
    [XLabl, YLabl, Titlename] = findRightLabel((propertyID), timeUnit, lengthUnit, propertynames);

     %% Plot the data
   minv = min(data(:,2));
   maxv = max(data(:,2));
   if autoBins
       binNum = round(sqrt(size(data,1)));
       edges = linspace(minv, maxv, binNum);
       newBinNmbrs = binNum;
   else
       edges = linspace(minv, maxv, binNumbers);
       newBinNmbrs = binNumbers;
   end
   his = histogram(Axes, data(:,2), edges);
   hixMaxValue = max(his.Values);
   xlim(Axes, [minv maxv]);
   title(Axes, Titlename);
   xlabel(Axes, YLabl);
   ylabel(Axes, "Counts");
   
   %% Decide if we fit or not
   if performFit
       %% perform fit
       pdGauss = fitdist(data(:,2), "Normal");
       pdKernel = fitdist(data(:,2), "Kernel", "Width", []);

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
       gaussCi95 = paramci(pdGauss);
       kernelCi95 = [0,0]; %not possible for kernel
       gaussNLL = negloglik(pdGauss);
       kernelNLL = negloglik(pdKernel);
       gaussDat = [median(pdGauss), mean(pdGauss), std(pdGauss), (std(pdGauss)/sqrt(size(data,1))), var(pdGauss), gaussNLL];
       kernelDat = [median(pdKernel), mean(pdKernel), std(pdKernel), (std(pdKernel)/sqrt(size(data,1))), var(pdKernel), kernelNLL];
   else
       gaussDat = [];
       kernelDat = [];
   end

end

%at the end so it is not so disruptive...
function [XLabl, YLabl, Titlename] = findRightLabel(ID, time, length, propertynames)
    %allpropertyNames = ["TrackID", "Timestep", "X-Pos", "Y-Pos", "Z-Pos", "Jumpdistance", "Ext-Segment-D", "Ext-Segment-DErr", "Ext-Segment-DR", 
    %"Ext-Segment-MeanJumpDistance", "Ext-Segment-MSD", "Ext-Segment-MSDErr", "Motiontype", "SegmentStart", "SegmentLifetime", "NumberOfSegmentsInTrack", 
    % DSPT-Alpha", "DSPT-D", "DSPT-Efficiency", "DSPT-logEfficiency", "DSPT-FractalDim", "DSPT-Gaussiantiy", "DSPT-Kurtosis", "DSPT-MSDratio", "DSPT-Trappedness", 
    % "DSPT-T0", "DSPT-T1", "DSPT-T2", "DSPT-T3", "DSPT-Lifetime", "DSPT-Length", "DSPT-AvgStepLength", "DSPT-AvgMSD", "DSPT-AvgDP", "DSPT-corrDP", "DSPT-signDP", 
    % "DSPT-SumSteplength", "DSPT-MinSteplength", "DSPT-MaxSteplength", "DSPT-BroadnessSteplength", "DSPT-Speed", "DSPT-Covariance", "DSPT-FractionSlow", 
    % "DSPT-FractionFast", "DSPT-Volume", "DSPT-Percent-NormalDiff", "DSPT-Percent-DirectedMot", "DSPT-Percent-ConfinedDiff", "DSPT-Percent-SuperDiff", 
    % "DSPT-NumberChangepoints", "DSPT-Instant-MSD-D", "DSPT-meanSequence", "DSPT-medianSequence", "maxSequence", "DSPT-minSequence", "DSPT-stdSequence", 
    % "DSPT-SimSeq"];
    %allpropertyIDS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,0,0,0,0,0,23,24,0,0,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
    switch ID
        case 1
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(1));
            Titlename = propertynames(1);
        case 2
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(2));
            Titlename = propertynames(2);
        case 3
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(3), length);
            Titlename = propertynames(3);
        case 4
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(4), length);
            Titlename = propertynames(4);
        case 5
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(5), length);
            Titlename = propertynames(5);
        case 6
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(6), length);
            Titlename = propertynames(6);
        case 7
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(7), "µm²/s");
            Titlename = propertynames(7);
        case 8
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(8), "µm²/s");
            Titlename = propertynames(8);
        case 9
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s ", propertynames(9));
            Titlename = propertynames(9);
        case 10
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(10), length);
            Titlename = propertynames(10);
        case 11
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(11), "µm²/s²");
            Titlename = propertynames(11);
        case 12
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(12), "µm²/s²");
            Titlename = propertynames(12);
        case 13
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(13));
            Titlename = propertynames(13);
        case 14
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(14));
            Titlename = propertynames(14);
        case 15
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(15), time);
            Titlename = propertynames(15);
        case 16
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(16), "steps");
            Titlename = propertynames(16);
        case 17
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(17));
            Titlename = propertynames(17);
        case 18
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(18), "a.u.");
            Titlename = propertynames(18);
        case 19
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(19), length);
            Titlename = propertynames(19);
        case 20
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(20));
            Titlename = propertynames(20);
        case 21
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(21));
            Titlename = propertynames(21);
        case 22
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(22));
            Titlename = propertynames(22);
        case 23
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(23));
            Titlename = propertynames(23);
        case 24
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(24), join([length "²/" time],"") );
            Titlename = propertynames(24);
        case 25
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(25));
            Titlename = propertynames(25);
        case 26
                                                XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(26));
            Titlename = propertynames(26);
        case 27
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(27));
            Titlename = propertynames(27);
        case 28
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(28));
            Titlename = propertynames(28);
        case 29
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(29));
            Titlename = propertynames(29);
        case 30
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(30));
            Titlename = propertynames(30);
        case 31
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(31));
            Titlename = propertynames(31);
        case 32
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(32));
            Titlename = propertynames(32);
        case 33
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(33), length);
            Titlename = propertynames(33);
        case 34
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(34));
            Titlename = propertynames(34);
        case 35
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(35));
            Titlename = propertynames(35);
        case 36
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(36));
            Titlename = propertynames(36);
        case 37
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(37));
            Titlename = propertynames(37);
        case 38
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(38), "steps");
            Titlename = propertynames(38);
        case 39
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(39), "steps");
            Titlename = propertynames(39);
        case 40
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(40), length);
            Titlename = propertynames(40);
        case 41
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(41), join([length "²/" time "²"],""));
            Titlename = propertynames(41);
        case 42
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(42));
            Titlename = propertynames(42);
        case 43
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(43));
            Titlename = propertynames(43);
        case 44
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(44));
            Titlename = propertynames(44);
        case 45
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(45), length);
            Titlename = propertynames(45);
        case 46
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(46), length);
            Titlename = propertynames(46);
        case 47
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(47), length);
            Titlename = propertynames(47);
        case 48
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(48), length);
            Titlename = propertynames(48);
        case 49
                        XLabl = sprintf("Time [%s]", time); %what is this
            YLabl = sprintf("%s [%s]", propertynames(49), join([length,"/",time],""));
            Titlename = propertynames(49);
        case 50
                        XLabl = sprintf("Time [%s]", time); %what is this
            YLabl = sprintf("%s", propertynames(50));
            Titlename = propertynames(50);
        case 51
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(51), "%");
            Titlename = propertynames(51);
        case 52
                        XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(52), "%");
            Titlename = propertynames(52);
        case 53
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(53), join([length, "³"],""));
            Titlename = propertynames(53);
        case 54
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(54), "%");
            Titlename = propertynames(54);
        case 55
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(55), "%");
            Titlename = propertynames(55);
        case 56
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(56), "%");
            Titlename = propertynames(56);
        case 57
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(57), "%");
            Titlename = propertynames(57);
        case 58
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(58));
            Titlename = propertynames(58);
        case 59
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s [%s]", propertynames(59), join([length "²/" time "²"],""));
            Titlename = propertynames(59);
        case 60
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(60));
            Titlename = propertynames(60);
        case 61
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(61));
            Titlename = propertynames(61);
        case 62
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(62));
            Titlename = propertynames(62);
        case 63
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(63));
            Titlename = propertynames(63);
        case 64
                                    XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(64));
            Titlename = propertynames(64);
        case 65
                                                XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(65));
            Titlename = propertynames(65);
            
    end
end