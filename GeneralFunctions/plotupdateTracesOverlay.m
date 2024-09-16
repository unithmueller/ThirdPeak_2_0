function [lineColors, spotColors, locPrecCb] = plotupdateTracesOverlay(ax1, Traces, flag3d, On, type, frame, scale, offsetx, offsety, offsetz, spotsize, lineC, spotC, zslicePos, zScalingValue, isImage3D, showLocPrec, LocPrecColorbar, unitUsed)
%function to plot the tracks on top of the micrscopy image. can select to
%either be in 2d or 3d tracks, as well as if the localisations or a
%connected track should be shown. 
%Will need to figure out how to make the units of the tracks an the pixels
%work so it fits...
    arguments
        ax1 %UI axes
        Traces (:,:) %trackdata
        flag3d (1,1) %if 3d or not
        On %if overlay plottet or not
        type %traces or spots plottet
        frame %current frame in slider
        scale (1,1)
        offsetx (1,1)
        offsety (1,1)
        offsetz (1,1)
        spotsize (1,1)
        lineC
        spotC
        zslicePos
        zScalingValue
        isImage3D
        showLocPrec
        LocPrecColorbar
        unitUsed
    end
    
    if isempty(LocPrecColorbar)
        locPrecCb = [];
    else
        locPrecCb = LocPrecColorbar;

    end

    if On == 0
        %nothing
        lineColors = lineC;
        spotColors = spotC;
    else
        %get frame relevant data first
        spots = Traces(Traces(:,2) == frame,:);
        if isempty(spots)
            lineColors = lineC;
            spotColors = spotC;
            return
        end
        trackIds = Traces(Traces(:,2) == frame,1);
        if isempty(trackIds)
            lineColors = lineC;
            spotColors = spotC;
            return
        end
        %for 3d stack data we filter for the current z slice
        if isImage3D
            spots = Traces((Traces(:,2) == frame) & (Traces(:,5)*zScalingValue == zslicePos),:);
            if isempty(spots)
                lineColors = lineC;
                spotColors = spotC;
                return
            end
            trackIds = Traces((Traces(:,2) == frame) & (Traces(:,5)*zScalingValue == zslicePos),1);
            if isempty(trackIds)
                lineColors = lineC;
                spotColors = spotC;
            return
            end
        end
        %need for loop...
        tracks = cell(size(trackIds,1),1);
        for i = 1:size(trackIds,1)
            tracks{i,1} = Traces(Traces(:,1) == trackIds(i),:);
        end
        %change scaling to match image data
        spots(:,3:4) = spots(:,3:4)*scale;
        spots(:,3) = spots(:,3)+offsetx;
        spots(:,4) = spots(:,4)+offsety;
        spots(:,5) = spots(:,5)+offsetz;
        
        for i = 1:size(trackIds,1)
            tracks{i,1}(:,3:4) = tracks{i,1}(:,3:4)*scale;
            tracks{i,1}(:,3) = tracks{i,1}(:,3)+offsetx;
            tracks{i,1}(:,4) = tracks{i,1}(:,4)+offsety;
            tracks{i,1}(:,5) = tracks{i,1}(:,5)+offsetz;
        end
        
        %need to plot
        hold(ax1, "on");
        handles = struct();
        %if isempty(spotC)
        %allocate some save space for the colors
        spotColors = cell(size(spotC,1) + size(trackIds,1),2);
        lineColors = cell(size(lineC,1) + size(trackIds,1),2);
        %add the old colors
        if ~isempty(spotC)
            spotColors(1:size(spotC,1),:) = spotC;
        end
        if ~isempty(lineC)
            lineColors(1:size(lineC,1),:) = lineC;
        end
        %carry a counter for new additions
        colorCounterSpot = size(spotC,1)+1;
        colorCounterLine = size(lineC,1)+1;

        %spotColors = spotC;
        %lineColors = lineC;
        PrecisionRange = max(Traces(:,6))-min(Traces(:,6));
        PrecisionArray = linspace(max(Traces(:,6)),min(Traces(:,6)),256);
        PrecisionColorMap = (colormap());

        for i = 1:size(trackIds)
            %for every track in the current frame
            switch type
                case 0 %only spots
                    if flag3d == 0 %2D
                        curTrackID = trackIds(i);
                        try
                            IDs = spotColors(:,1);
                            IDs = cell2mat(IDs);
                        catch
                            IDs = [-1; -1];
                        end
                        if showLocPrec
                            PrecColorPos = spots(spots(:,1) == trackIds(i),6);
                            PrecColors = [];
                            for j = 1:size(PrecColorPos,1)
                                [~,PrecIdx] = min(abs(PrecisionArray-PrecColorPos(j)));
                                PrecColor = PrecisionColorMap(PrecIdx,:);
                                PrecColors(j,:) = PrecColor;
                            end
                            
                            h = scatter(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spotsize, PrecColors, 'filled');
                            if isempty(LocPrecColorbar)
                                cbFig = figure();
                                cbAx = axes(cbFig);
                                h1 = scatter(cbAx, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spotsize, PrecColors, 'filled');
                                colormap(cbAx, flip(PrecisionColorMap));
                                clim(cbAx, [min(PrecisionArray) max(PrecisionArray)]);
                                colbar = colorbar(cbAx);
                                colbar.Label.String = sprintf("Localization precision [%s]", unitUsed);
                                axChildren = cbAx.Children;
                                title(cbAx, "Localization precision");
                                set(axChildren, "Visible", "off");
                                cbAx.Visible = "off";
                                LocPrecColorbar = cbFig;
                                locPrecCb = LocPrecColorbar;
                            end
                        else
                            locPrecCb = [];
                            if any(IDs == curTrackID)
                                curColor = IDs == curTrackID;
                                curColor = spotColors{curColor,2};
                                h = scatter(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spotsize, curColor);
                            else
                                h = scatter(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spotsize);
                                spotColors{colorCounterSpot,1} = trackIds(i);
                                spotColors{colorCounterSpot,2} = h.CData;
                                colorCounterSpot = colorCounterSpot+1;
                            end
                        end
                    elseif flag3d == 1 %3D
                        curTrackID = trackIds(i);
                        try
                            IDs = spotColors(:,1);
                            IDs = cell2mat(IDs);
                        catch
                            IDs = [-1; -1];
                        end
                        if showLocPrec
                            PrecColorPos = spots(spots(:,1) == trackIds(i),6);
                            PrecColors = [];
                            for j = 1:size(PrecColorPos,1)
                                [~,PrecIdx] = min(abs(PrecisionArray-PrecColorPos(j)));
                                PrecColor = PrecisionColorMap(PrecIdx,:);
                                PrecColors(j,:) = PrecColor;
                            end
                            h = scatter3(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spots(spots(:,1) == trackIds(i),5), spotsize, PrecColors, 'filled');
                            if isempty(LocPrecColorbar)
                                cbFig = figure();
                                cbAx = axes(cbFig);
                                h1 = scatter3(cbAx, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spots(spots(:,1) == trackIds(i),5), spotsize, PrecColors, 'filled');
                                colormap(cbAx, flip(PrecisionColorMap));
                                clim(cbAx, [min(PrecisionArray) max(PrecisionArray)]);
                                colorbar(cbAx);
                                colbar.Label.String = sprintf("Localization precision [%s]", unitUsed);
                                axChildren = cbAx.Children;
                                title(cbAx, "Localization precision");
                                axChildren = cbAx.Children;
                                set(axChildren, "Visible", "off");
                                cbAx.Visible = "off";
                                LocPrecColorbar = cbFig;
                                locPrecCb = LocPrecColorbar;
                            end
                        else
                            locPrecCb = [];
                            if any(IDs == curTrackID)
                                curColor = IDs == curTrackID;
                                curColor = spotColors{curColor,2};
                                h = scatter3(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spots(spots(:,1) == trackIds(i),5), spotsize, curColor);
                            else
                                h = scatter3(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spots(spots(:,1) == trackIds(i),5), spotsize);
                                spotColors{colorCounterSpot,1} = trackIds(i);
                                spotColors{colorCounterSpot,2} = h.CData;
                                colorCounterSpot = colorCounterSpot+1;
                            end
                        end
                    end

                case 1 %tracks and spots
                    %for every Track
                    curTrack = tracks{i,1};
                        if flag3d == 0
                            curTrackID = trackIds(i);
                            try
                                IDs = lineColors(:,1);
                                IDs = cell2mat(IDs);
                            catch
                                IDs = [-1; -1];
                            end
                            if any(IDs == curTrackID)
                                curColor = IDs == curTrackID;
                                curColor = lineColors{curColor,2};
                                h = plot(ax1, curTrack(:,3), curTrack(:,4), 'Color',curColor);
                                handles.line(i) = h;
                                color = handles.line(i).Color;
                                scatter(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spotsize, color, 'filled');
                            else
                                h = plot(ax1, curTrack(:,3), curTrack(:,4));
                                handles.line(i) = h;
                                color = handles.line(i).Color;
                                lineColors{colorCounterLine,1} = trackIds(i);
                                lineColors{colorCounterLine,2} = color;
                                colorCounterLine = colorCounterLine+1;
                                scatter(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spotsize, color, 'filled');
                            end
                        elseif flag3d == 1
                            curTrackID = trackIds(i);
                            try
                                IDs = lineColors(:,1);
                                IDs = cell2mat(IDs);
                            catch
                                IDs = [-1; -1];
                            end
                            if any(IDs == curTrackID)
                                curColor = IDs == curTrackID;
                                curColor = lineColors{curColor,2};
                               h = plot3(ax1, curTrack(:,3), curTrack(:,4), curTrack(:,5), 'Color', curColor);
                                handles.line(i) = h;
                                color = handles.line(i).Color;
                                scatter3(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spots(spots(:,1) == trackIds(i),5), spotsize, color, 'filled');
                            else
                                h = plot3(ax1, curTrack(:,3), curTrack(:,4), curTrack(:,5));
                                handles.line(i) = h;
                                color = handles.line(i).Color;
                                lineColors{colorCounterLine,1} = trackIds(i);
                                lineColors{colorCounterLine,2} = color;
                                colorCounterLine = colorCounterLine+1;
                                scatter3(ax1, spots(spots(:,1) == trackIds(i),3), spots(spots(:,1) == trackIds(i),4), spots(spots(:,1) == trackIds(i),5), spotsize, color, 'filled');
                            end
                        end
             end
        end
        %remove empty color placeholder
        spotColors = spotColors(~cellfun(@isempty, spotColors(:,1)),:);
        lineColors = lineColors(~cellfun(@isempty, lineColors(:,1)),:);
    end
end