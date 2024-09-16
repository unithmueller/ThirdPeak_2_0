function makeScalarMapFromProperty(app, data, flag3d, superpxsize, NmbrTimeWins, timeUnit, lengthUnit, superZsize, unit, propertyID, propertyNames, allPropertyIDs)
%generates a supermap of velocity data in 2d or 3d using the given pxsize
    %need to convert data always into unit/frame
    type = split(convertCharsToStrings(app.ImportSettingsStruct.customUnits.DataType),"/");
            comparison = strcmp(type, ["Pixel"; "Frame"]);
            if comparison %pixel/frame
                %if pixel/frame, covert pixel to unit
                data(:,3:5) = data(:,3:5)*app.ImportSettingsStruct.customUnits.Pixelsize;
            elseif comparison(1) & (~comparison(2))%Pixel/Unit
                %if pixel/unit, must convert unit to frame and pixel to
                %unit
                time = data(:,2);
                time = round(time/app.ImportSettingsStruct.customUnits.Timestep,0);
                data(:,2) = time;
                data(:,3:5) = data(:,3:5)*app.ImportSettingsStruct.customUnits.Pixelsize;
            elseif (~comparison(1)) & comparison(2)%Unit/Frame
                %if unit/frame, do nothing
            else
                %if unit/unit, must convert time to frame
                time = data(:,2);
                time = round(time/app.ImportSettingsStruct.customUnits.Timestep,0);
                data(:,2) = time;
            end
    
    %Set the boundaries of the scalar map based on the spatial properties
    minX = min(data(:,3));
    maxX = max(data(:,3));
    minY = min(data(:,4));
    maxY = max(data(:,4));
    minZ = min(data(:,5));
    maxZ = max(data(:,5));
    TimeWinDur = floor(max(data(:,2))/NmbrTimeWins);
    dx = superpxsize;
    dz = superZsize;
    
    tSegData = cell(1,NmbrTimeWins);
    %Segment the data in time if necessary
    if NmbrTimeWins == 1
        %single win
        tSegData{1,1} = data;
    else
        %many win
        begin = find(min(data(:,2)));
        for i = 1:NmbrTimeWins
            ending = begin+TimeWinDur;
            tSegData{1,i} = data(begin:ending,:);
            begin = ending+1;
        end
    end
    
    %set the edge values
    xEdgeNumbers = ceil((maxX-minX)/dx);
    xEdge = minX + dx * [0:1:xEdgeNumbers]; 
    yEdgeNumbers = ceil((maxY-minY)/dx);
    yEdge = minY + dx * [0:1:yEdgeNumbers];

    %make cell array to store the data
    if flag3d
        map = cell(ceil((maxX-minX)/dx), ceil((maxY-minY)/dx), ceil((maxZ-minZ)/dz));
        zEdgeNumbers = ceil((maxZ-minZ)/dz);
        zEdge = minZ + dz*[0:1:zEdgeNumbers];
    else
        map = cell(ceil((maxX-minX)/dx), ceil((maxY-minY)/dx));
    end

    maps = cell(1,NmbrTimeWins);

    %seperate into bins
    if flag3d
        for j = 1:NmbrTimeWins
            tmpMap = map;
            [~,~,~,binX,binY] = histcounts2(tSegData{1,j}(:,3), tSegData{1,j}(:,4), xEdge, yEdge); 
            [~,~,binZ] = histcounts(tSegData{1,j}(:,5),zEdge);
            %make the cell array with the sptially sorted data
            for i = 1:size(binX,1)
                if ((propertyID > 53) & (propertyID < 58))
                    [~,tmpdat] = max(tSegData{1,j}(i,54:57));
                    tmpMap{binX(i), binY(i), binZ(i)} = [tmpMap{binX(i), binY(i)}; tmpdat];
                else
                    tmpMap{binX(i), binY(i), binZ(i)} = [tmpMap{binX(i), binY(i), binZ(i)}; tSegData{1,j}(i,propertyID)];
                end
            end
            if propertyID == 14 | ((propertyID > 53) & (propertyID < 58))
                tmpMap = cellfun(@mode, tmpMap, "UniformOutput",true);
            else
                tmpMap = cellfun(@mean, tmpMap, "UniformOutput",true);
            end
            maps{1,j} = tmpMap;
        end
    else
        for j = 1:NmbrTimeWins
            tmpMap = map;
            [~,~,~,binX,binY] = histcounts2(tSegData{1,j}(:,3), tSegData{1,j}(:,4), xEdge, yEdge);
             %make the cell array with the sptially sorted data
            for i = 1:size(binX,1)
                if ((propertyID > 53) & (propertyID < 58))
                    [~,tmpdat] = max(tSegData{1,j}(i,54:57));
                    tmpMap{binX(i), binY(i)} = [tmpMap{binX(i), binY(i)}; tmpdat];
                else
                    tmpMap{binX(i), binY(i)} = [tmpMap{binX(i), binY(i)}; tSegData{1,j}(i,propertyID)];
                end
            end
            %might need other functions but could switch them out if needed
            if propertyID == 14 | ((propertyID > 53) & (propertyID < 58))
                tmpMap = cellfun(@mode, tmpMap, "UniformOutput",true);
            else
                tmpMap = cellfun(@mean, tmpMap, "UniformOutput",true);
            end
            maps{1,j} = tmpMap;
        end
    end

    a = 1;
    %% now the plotting
    %Prepare the labels
    propName = propertyNames(propertyID);
    [XLabl, YLabl, Titlename] = findRightLabel(propertyID, timeUnit, lengthUnit, propertyNames);

    %Do the plotting
    for i = 1:NmbrTimeWins
        data = flip(rot90(maps{1,i}));
        %hold off
        figure('Name', sprintf('Property Map - %s - TimeWindow %d', Titlename, i));
        if flag3d
            %3d
            [X, Y] = meshgrid(1:size(data,2),1:size(data,1));
            for Z = 1:size(data,3)
                hold on;
                C(Y,X) = data(Y,X,Z);
                C = C;
                zpos = zeros(size(data,1),size(data,2));
                zpos(:) = Z;
                s = surface(X,Y,zpos,C);
                s.EdgeColor = 'none';
                alpha(s,"color");
                alpha(s,"scaled");
                %hold on
            end
            c = colorbar;
            c.Label.String = sprintf('%s',YLabl);
            xlabel("X Position [px]");
            ylabel("Y Position [px]");
            zlabel("Z Position [px]");
            title(sprintf("%s%d" , Titlename, i));
            hold(gca,"off");
            clear C;
        else
            %2d
            %figure('Name','Velocity Map');

            data = data;
            [X, Y] = meshgrid(1:size(data,2),1:size(data,1));
            Z = zeros(size(data));
            s = surface(X,Y,Z,data(:,:));
            s.EdgeColor = 'none';
            %alpha(s,"color");
            %alpha(s,"scaled");
            c = colorbar;
            c.Label.String = sprintf('%s',YLabl);
            xlabel("X Position [px]");
            ylabel("Y Position [px]");
            title(sprintf("%s%d" , Titlename, i));
            hold off
        end
        if propertyID == 14
            custMap =[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
            colormap(custMap);
            clim([1 5]);
            c.Ticks = ([1 2 3 4 5]);
            c.TickLabels = {'none','immobile','diffusion','directed','directedDif'};
        elseif (propertyID >53) & (propertyID < 58)
            custMap =[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
            colormap(custMap);
            clim([1 5]);
            c.Ticks = ([1 2 3 4 5]);
            c.TickLabels = {'diffusion','directed','confined','anomal','none'};
        end
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
            YLabl = "";%sprintf("%s [%s]", propertynames(54), "%");
            Titlename = propertynames(54);
        case 55
            XLabl = sprintf("Time [%s]", time);
            YLabl = "";%sprintf("%s [%s]", propertynames(55), "%");
            Titlename = propertynames(55);
        case 56
            XLabl = sprintf("Time [%s]", time);
            YLabl = "";%sprintf("%s [%s]", propertynames(56), "%");
            Titlename = propertynames(56);
        case 57
            XLabl = sprintf("Time [%s]", time);
            YLabl = "";%sprintf("%s [%s]", propertynames(57), "%");
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
        case 66
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(66));
            Titlename = propertynames(66);
        case 67
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(67));
            Titlename = propertynames(67);
        case 68
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(68));
            Titlename = propertynames(68);
        case 69
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(69));
            Titlename = propertynames(69);
        case 70
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(70));
            Titlename = propertynames(70);
        case 71
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(71));
            Titlename = propertynames(71);
        case 72
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(72));
            Titlename = propertynames(72);
        case 73
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(73));
            Titlename = propertynames(73);
        case 74
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(74));
            Titlename = propertynames(74);
        case 75
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(75));
            Titlename = propertynames(75);
        case 76
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(76));
            Titlename = propertynames(76);
        case 77
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(77));
            Titlename = propertynames(77);
        case 78
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(78));
            Titlename = propertynames(78);
        case 79
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(79));
            Titlename = propertynames(79);
        case 80
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(80));
            Titlename = propertynames(80);
        case 81
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(81));
            Titlename = propertynames(81);
        case 82
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(82));
            Titlename = propertynames(82);
        case 83
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(83));
            Titlename = propertynames(83);
        case 84
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(84));
            Titlename = propertynames(84);
        case 85
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(85));
            Titlename = propertynames(85);
        case 86
            XLabl = sprintf("Time [%s]", time);
            YLabl = "°";
            Titlename = propertynames(86);
        case 87
            XLabl = sprintf("Time [%s]", time);
            YLabl = "°";
            Titlename = propertynames(87);
        case 88
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(88));
            Titlename = propertynames(88);
        case 89
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(89));
            Titlename = propertynames(89);
        case 90
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(90));
            Titlename = propertynames(90);
        case 91
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(91));
            Titlename = propertynames(91);
        case 92
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(92));
            Titlename = propertynames(92);
        case 93
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(93));
            Titlename = propertynames(93);
        case 94
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(94));
            Titlename = propertynames(94);
        case 95
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(95));
            Titlename = propertynames(95);
        case 96
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(96));
            Titlename = propertynames(96);
        case 97
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(97));
            Titlename = propertynames(97);
        case 98
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(98));
            Titlename = propertynames(98);
        case 99
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(99));
            Titlename = propertynames(99);
        case 100
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(100));
            Titlename = propertynames(100);
        case 101
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(101));
            Titlename = propertynames(101);
        case 102
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(102));
            Titlename = propertynames(102);
        case 103
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(103));
            Titlename = propertynames(103);
        case 104
            XLabl = sprintf("Time [%s]", time);
            YLabl = sprintf("%s", propertynames(104));
            Titlename = propertynames(104);
    end
            
end
