function makeDiffusionMap(app, data, flag3d, superpxsize, NmbrTimeWins, minAcqDT, nump, filter, filtersize, superZsize, unit, expTime, locPrecXY, locPrecZ)
%generates a supermap of diffusion data in 2d or 3d using the given pxsize
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

    trajensemb = TrajectoryEnsemble(data(:,1:5), minAcqDT, flag3d);

    maxXY = max(data(:,3:4));
    minZ = min(data(:,5));
    maxZ = max(data(:,5));

    TimeWinDur = ceil(max(data(:,2))/NmbrTimeWins);
    
    tws = TimeWindows.empty;
    if NmbrTimeWins == 1
        %single win
        tmptws = TimeWindows.generateSingleWindow(trajensemb);
        tws(end+1) = tmptws;
    else
        %many win
        begin = min(data(:,2));
        for i = 1:NmbrTimeWins
            ending = begin+TimeWinDur;
            tmpwin = TimeWindows(TimeWinDur, 0, begin, ending);
            begin = ending+1;
            tws(end+1) = tmpwin;
        end
    end
    
    trajensembwins = TrajectoryEnsembleWindows(trajensemb, tws);
    
    params = DiffusionParameters(superpxsize, nump, filter, filtersize, flag3d, superZsize, expTime, locPrecXY, locPrecZ, minAcqDT);
    
    scalmapwins = ScalarMapWindows.gen_realdiffusion_maps(trajensembwins, params);
    a = 1;
    units = split(unit, "/");
    %now the plotting
    for i = 1:size(scalmapwins.windows,2)
        data = flip(rot90(scalmapwins.windows(1,i).data));
        hold off
        figure('Name', sprintf('Diffusion Map - TimeWindow %d', i));
        if flag3d
            %3d
            [X, Y] = meshgrid(1:size(data,2),1:size(data,1));
            for Z = 1:size(data,3)
                hold on;
                C(Y,X) = data(Y,X,Z);
                %C = C/(6*minAcqDT);
                C = C/(6);
                zpos = zeros(size(data,1),size(data,2));
                zpos(:) = Z;
                s = surface(X,Y,zpos,C);
                s.EdgeColor = 'none';
                alpha(s,"color");
                alpha(s,"none");
                %hold on
            end
            c = colorbar;
            txt = join([units(1) "²/" units(2)],"");
            c.Label.String = sprintf('Diffusion coefficient [%s]',txt);
            xlabel("X Position [px]");
            ylabel("Y Position [px]");
            zlabel("Z Position [px]");
            hold(gca,"off");
            clear C;
        else
            %2d
            %figure('Name','Velocity Map');

            %data = data/(4*minAcqDT); %from nm/timestep to nm/timeunitstep
            data = data/(4); %from nm/timestep to nm/timeunitstep
            [X, Y] = meshgrid(1:size(data,2),1:size(data,1));
            Z = zeros(size(data));
            s = surface(X,Y,Z,data(:,:));
            s.EdgeColor = 'none';
            %alpha(s,"color");
            %alpha(s,"scaled");
            c = colorbar;
            txt = join([units(1) "²/" units(2)],"");
            c.Label.String = sprintf('Diffusion coefficient [%s]',txt);
            xlabel("X Position [px]");
            ylabel("Y Position [px]");
            hold off
        end
    end     
end
