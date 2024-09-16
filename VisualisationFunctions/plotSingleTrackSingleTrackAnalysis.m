function plotSingleTrackSingleTrackAnalysis(Axes, Track, is3d, timepoint, lengthUnit, isPixel)
    %Function to plot the selected single track, use slider to show
    %interactively plot it. Use Color to classify the diffusion
    %properties/segments from DeepSPT or Swift if available.
    PlotLineWidth = 2;
    ScatterPointSize = 35;

    %Determine Data availability, 2 for Deep, 1 for Swft, 0 for nothing
    if size(Track,2) > 23
        DataAvailable = 2;
    elseif size(Track,2) > 5
        if mean(Track(:,6:22),"all") ~= 0
            DataAvailable = 1;
        else
            DataAvailable = 0;
        end
    else
        DataAvailable = 0;
    end

    %Determine dimension
    if is3d == "3D"
        is3d = 1;
    else
        is3d = 0;
    end

    %Determine relevant data
    TrackBasic = Track(:,2:5);
    if DataAvailable == 2
        %DeepSPT
        [maxVals, maxIDX] = max(Track(:,54:57),[],2);
        % percNorm = Track(:,46);
        % percDir = Track(:,47);
        % percConfin = Track(:,48);
        % percAnom = Track(:,49);
    elseif DataAvailable == 1
        motiontype = Track(:,14);
        translationArray = [5, 3, 1, 2, 4];
        % none = 1; immob = 2, diff = 3, direct = 4, dirdiff = 5
    end
    
    %make a clean slate
    cla(Axes,'reset');
    %Plot the data
    % 1=normdiff, 2=Direct, 3=Confined, 4=Anomal, 5=none
    customColormap = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
    if is3d
        if DataAvailable == 2
            initIDX = 1;
            hold(Axes, "on");
            for stps = 2:timepoint
                tmpPoints = TrackBasic(initIDX:stps,:);
                tmpType = maxIDX(initIDX);
                tmpPrec = maxVals(initIDX);
                tmpColor = customColormap(tmpType,:);
                tmpColor(end+1) = tmpPrec;
                p = plot3(Axes, tmpPoints(:,2), tmpPoints(:,3), tmpPoints(:,4), "Color", tmpColor, "LineWidth", PlotLineWidth);
                initIDX = initIDX+1;
                if stps==timepoint
                    tmpColor = customColormap(tmpType,:);
                    scatter3(Axes, tmpPoints(end,2), tmpPoints(end,3), tmpPoints(end,4), ScatterPointSize, tmpColor, "filled");
                end
            end
            hold(Axes, "off");
        elseif DataAvailable == 1
            initIDX = 1;
            hold(Axes, "on");
            for stps = 2:timepoint
                tmpPoints = TrackBasic(initIDX:stps,:);
                tmpType = motiontype(initIDX);
                tmpType = translationArray(tmpType);
                plot3(Axes, tmpPoints(:,2), tmpPoints(:,3), tmpPoints(:,4), "Color", customColormap(tmpType,:), "LineWidth", PlotLineWidth);
                initIDX = initIDX+1;
                if stps==timepoint
                    scatter3(Axes, tmpPoints(end,2), tmpPoints(end,3), tmpPoints(end,4), ScatterPointSize, customColormap(tmpType,:), "filled");
                end
            end
            hold(Axes, "off");
        else
            plot3(Axes, TrackBasic(1:timepoint,2), TrackBasic(1:timepoint,3), TrackBasic(1:timepoint,4), "Color", customColormap(5,:), "LineWidth", PlotLineWidth);
            hold(Axes, "on");
            scatter3(Axes, TrackBasic(timepoint,2), TrackBasic(timepoint,3), TrackBasic(timepoint,4), ScatterPointSize, customColormap(5,:), "filled");
            hold(Axes, "off");
        end
    else
        if DataAvailable == 2
            initIDX = 1;
            hold(Axes, "on");
            for stps = 2:timepoint
                tmpPoints = TrackBasic(initIDX:stps,:);
                tmpType = maxIDX(initIDX);
                tmpPrec = maxVals(initIDX);
                tmpColor = customColormap(tmpType,:);
                tmpColor(end+1) = tmpPrec;
                p = plot(Axes, tmpPoints(:,2), tmpPoints(:,3), "Color", tmpColor, "LineWidth", PlotLineWidth);
                initIDX = initIDX+1;
                if stps==timepoint
                    tmpColor = customColormap(tmpType,:);
                    scatter(Axes, tmpPoints(end,2), tmpPoints(end,3), ScatterPointSize, tmpColor, "filled");
                end
            end
            hold(Axes, "off");
        elseif DataAvailable == 1
            initIDX = 1;
            hold(Axes, "on");
            for stps = 2:timepoint
                tmpPoints = TrackBasic(initIDX:stps,:);
                tmpType = motiontype(initIDX);
                tmpType = translationArray(tmpType);
                plot(Axes, tmpPoints(:,2), tmpPoints(:,3), "Color", customColormap(tmpType,:), "LineWidth", PlotLineWidth);
                initIDX = initIDX+1;
                if stps==timepoint
                    scatter(Axes, tmpPoints(end,2), tmpPoints(end,3), ScatterPointSize, customColormap(tmpType,:), "filled");
                end
            end
            hold(Axes, "off");
        else
            plot(Axes, TrackBasic(1:timepoint,2), TrackBasic(1:timepoint,3), "Color", customColormap(5,:), "LineWidth", PlotLineWidth);
            hold(Axes, "on");
            scatter(Axes, TrackBasic(timepoint,2), TrackBasic(timepoint,3), ScatterPointSize, customColormap(5,:), "filled");
            hold(Axes, "off");
        end
    end
    if isPixel
        xlabel(Axes, "X [px]");
        ylabel(Axes, "Y [px]");
    else
        xlabel(Axes, sprintf("X [%s]", lengthUnit));
        ylabel(Axes, sprintf("Y [%s]", lengthUnit));
    end
    title(Axes, sprintf("TrackID: %d", Track(1:1)))
    colormap(Axes, customColormap)
    colorbar(Axes, 'Ticks',[0.1,0.3,0.5,0.7,0.9], 'TickLabels', {'Normal', 'Directed', 'Confined', 'Anormal', 'None'});
    axis(Axes, "equal");
end