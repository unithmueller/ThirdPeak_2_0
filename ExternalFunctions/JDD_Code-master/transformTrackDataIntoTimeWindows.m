function [jumpDistances, x, y, z] = transformTrackDataIntoTimeWindows(TrackData, dimension, windowSize, isSlidingWindow, trackLengths,isPixel,timestep)
%Function to transform the data from within the software into the jump
%distances used for the JDD analyis. Handles sliding window and holes in
%the tracks
if isPixel
    factor = 1;
else
    factor = timestep;
end
if isSlidingWindow
    %make more out of longer tracks
    windowedTrackData = {};
    for i = 1:size(trackLengths,1)
        tmpTrack = TrackData(TrackData(:,1) == trackLengths(i), 2:5);
        tmpSteps = zeros((windowSize+1),3,(size(tmpTrack,1)-(windowSize)));
        for p = 1:(size(tmpTrack,1)-(windowSize))
            if tmpTrack(p+windowSize,1)-tmpTrack(p,1) == windowSize*factor
                %no hole, perfect step
                tmpSteps(:,:,p) = tmpTrack(p:p+windowSize,2:4);
            else
                %hole in step, need to solve the issue, try to find the
                %right step if possible
                for k = windowSize-1:-1:1
                    if tmpTrack(p+k,1)-tmpTrack(p,1) == windowSize*factor
                         tmpdat = tmpTrack(1:1+k+p,1:4);
                         holesPos = [];
                         for j = 1:size(tmpdat,1)-1
                             %determine holes
                             if tmpdat(j+1,1)-tmpdat(j,1) ~= factor
                                 holeSz = ((tmpdat(j+1,1)-tmpdat(j,1))/factor)-1;
                                 holesPos(j,1:2) = [j, holeSz];
                             end
                         end
                         holesPos(holesPos(:,2)==0,:) = [];
                         %generate placeholders to fill the void and
                         %emtpyness inside
                         placeholders = {};
                         for j = 1:size(holesPos,1)
                             placesToFill = holesPos(j,2);
                             xplaceH = linspace(tmpdat(holesPos(j,1),2),tmpdat(holesPos(j,1)+1,2),placesToFill+2);
                             xplaceH = xplaceH(2:end-1).';
                             yplaceH = linspace(tmpdat(holesPos(j,1),3),tmpdat(holesPos(j,1)+1,3),placesToFill+2);
                             yplaceH = yplaceH(2:end-1).';
                             zplaceH = linspace(tmpdat(holesPos(j,1),4),tmpdat(holesPos(j,1)+1,4),placesToFill+2);
                             zplaceH = zplaceH(2:end-1).';
                             placeH = [xplaceH, yplaceH, zplaceH];
                             placeholders{j} = placeH;
                         end

                         %put it together
                         putIndxs = [];
                         for j = 1:size(holesPos,1)
                             if j == 1
                                 putIndxs(j) = holesPos(j)+1;
                             else
                                 putIndxs(j) = holesPos(j)+1+sum(holesPos(1:j-1,2));
                             end
                         end

                         tmpdat = tmpdat(:,2:4);
                         for j = 1:size(putIndxs,2)
                              tmpdat = [tmpdat(1:putIndxs(j)-1,:); placeholders{j}; tmpdat(putIndxs(j):end,:)];
                         end
                         tmpSteps(:,:,j) = tmpdat(1:1+windowSize,:);
                    end
                end
            end
        end
        %cleanup leftovers
        tmpSteps(:,:,mean(mean(tmpSteps(:,:,1:end),1),2) == 0) = [];
        windowedTrackData{i} = tmpSteps;
    end
    totalTrackNumber = sum(cellfun(@(x) size(x,3), windowedTrackData));
    xdat = cellfun(@(x) reshape(x(:,1,:),windowSize+1,[]),windowedTrackData,'UniformOutput',false);
    x = cat(2, xdat{:});
    ydat = cellfun(@(x) reshape(x(:,2,:),windowSize+1,[]),windowedTrackData,'UniformOutput',false);
    y = cat(2, ydat{:});
    zdat = cellfun(@(x) reshape(x(:,3,:),windowSize+1,[]),windowedTrackData,'UniformOutput',false);
    z = cat(2, zdat{:});

    %Provide savelocation
    jumpDistances = zeros(size(totalTrackNumber,1),1);
    %Select the dimension
    if size(dimension,2) == 1
        switch dimension
            case 3
                dimVal = (x);
            case 4
                dimVal = (y);
            case 5
                dimVal = (z);
        end
        jumpDistances = sqrt((dimVal(end,:)-dimVal(1,:)).^2).';
    elseif size(dimension,2) == 2
        jumpDistances = sqrt((x(end,:)-x(1,:)).^2+(y(end,:)-y(1,:)).^2).';
    elseif size(dimension,2) == 3
        jumpDistances = sqrt((x(end,:)-x(1,:)).^2+(y(end,:)-y(1,:)).^2+(z(end,:)-z(1,:)).^2).';
    end

else
    %Only take the fixed window size
    %Provide savelocation
    jumpDistances = zeros(size(trackLengths,1),1);
    tmpSteps = zeros((windowSize+1),3,size(trackLengths,1));
    for i = 1:size(trackLengths,1)
        tmpTrack = TrackData(TrackData(:,1) == trackLengths(i), 2:5);
        if tmpTrack(1+windowSize,1)-tmpTrack(1,1) == windowSize*factor
            %no hole, perfect step
            tmpSteps(:,:,i) = tmpTrack(1:1+windowSize,2:4);
        else
            %hole in step, need to solve the issue, try to find the
            %right step if possible
            for k = windowSize-1:-1:1
                if tmpTrack(1+k,1)-tmpTrack(1,1) == windowSize*factor
                    tmpdat = tmpTrack(1:1+k,1:4);
                    holesPos = [];
                    for j = 1:(size(tmpdat,1)-1)
                        %determine holes
                        if tmpdat(j+1,1)-tmpdat(j,1) ~= factor
                            holeSz = ((tmpdat(j+1,1)-tmpdat(j,1))/factor)-1;
                            holesPos(j,1:2) = [j, holeSz];
                        end
                    end
                    holesPos(holesPos(:,2)==0,:) = [];
                    %generate placeholders to fill the void and
                    %emtpyness inside
                    placeholders = {};
                    for j = 1:size(holesPos,1)
                        placesToFill = holesPos(j,2);
                        xplaceH = linspace(tmpdat(holesPos(j,1),2),tmpdat(holesPos(j,1)+1,2),placesToFill+2);
                        xplaceH = xplaceH(2:end-1).';
                        yplaceH = linspace(tmpdat(holesPos(j,1),3),tmpdat(holesPos(j,1)+1,3),placesToFill+2);
                        yplaceH = yplaceH(2:end-1).';
                        zplaceH = linspace(tmpdat(holesPos(j,1),4),tmpdat(holesPos(j,1)+1,4),placesToFill+2);
                        zplaceH = zplaceH(2:end-1).';
                        placeH = [xplaceH, yplaceH, zplaceH];
                        placeholders{j} = placeH;
                    end

                    %put it together
                    putIndxs = [];
                    for j = 1:size(holesPos,1)
                        if j == 1
                            putIndxs(j) = holesPos(j)+1;
                        else
                            putIndxs(j) = holesPos(j)+1+sum(holesPos(1:j-1,2));
                        end
                    end
                    
                    tmpdat = tmpdat(:,2:4);
                    for j = 1:size(putIndxs,2)
                        tmpdat = [tmpdat(1:putIndxs(j)-1,:); placeholders{j}; tmpdat(putIndxs(j):end,:)];
                    end

                    tmpSteps(:,:,j) = tmpdat(1:1+windowSize,:);
                    
                end
            end
        end
    end

    totalTrackNumber = size(trackLengths,1);
    x = reshape(tmpSteps(:,1,:),windowSize+1,[]);
    y = reshape(tmpSteps(:,2,:),windowSize+1,[]);
    z = reshape(tmpSteps(:,3,:),windowSize+1,[]);

    %Provide savelocation
    jumpDistances = zeros(size(totalTrackNumber,1),1);
    %Select the dimension
    if size(dimension,2) == 1
        switch dimension
            case 3
                dimVal = (x);
            case 4
                dimVal = (y);
            case 5
                dimVal = (z);
        end
        jumpDistances = sqrt((dimVal(end,:)-dimVal(1,:)).^2).';
    elseif size(dimension,2) == 2
        jumpDistances = sqrt((x(end,:)-x(1,:)).^2+(y(end,:)-y(1,:)).^2).';
    elseif size(dimension,2) == 3
        jumpDistances = sqrt((x(end,:)-x(1,:)).^2+(y(end,:)-y(1,:)).^2+(z(end,:)-z(1,:)).^2).';
    end

end

