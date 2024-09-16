function TrackData = addSegmentToTracks(TrackData)
    %adds an additional ID at the end of the datasets to identify and
    %filter by track segments, not only whole tracks.
    datSize = size(TrackData{1,1});

    if datSize(2) == 5
        %only raw tracks, we dont know anything about segments
        for i = 1:size(TrackData,1)
            tmpdat = TrackData{i,1};
            tmpdat(:,end+1) = tmpdat(:,1);
            TrackData{i,1} = tmpdat;
        end
    elseif datSize(2) == 22
        %swift tracking, 22 = number of seg, 21 = start of seg with 1, stop
        %with 0
        for i = 1:size(TrackData,1)
            tmpdat = TrackData{i,1};
            newTmpDat = zeros(size(tmpdat,1), size(tmpdat,2)+1);
            bin = unique(tmpdat(:,1));
            counter = 1;
            for k = 1:size(bin,1)
                track = tmpdat(tmpdat(:,1) == bin(k),:);
                if mean(track(:,22) == 1)
                    track(:,end+1) = counter;
                    newTmpDat(tmpdat(:,1) == track(1,1),:) = track;
                    counter = counter+1;
                else
                    %more segments
                    trackSegsEndIdx = find(track(:,20) == 0);
                    segNumbers = zeros(size(track,1),1);
                    strt = 1;
                    for j = 1:length(trackSegsEndIdx)
                        segNumbers(strt:trackSegsEndIdx(j)) = counter;
                        strt = trackSegsEndIdx(j)+1;
                        counter = counter+1;
                    end
                    segNumbers(segNumbers == 0) = counter;
                    counter = counter +1;
                    track(:,end+1) = segNumbers;
                    newTmpDat(tmpdat(:,1) == track(1,1),:) = track;
                end
            end
            TrackData{i,1} = newTmpDat;
        end
    elseif datSize(2) == 65
        %deepspt
        for i = 1:size(TrackData,1)
            tmpdat = TrackData{i,1};
            newTmpDat = zeros(size(tmpdat,1), size(tmpdat,2)+1);
            bin = unique(tmpdat(:,1));
            counter = 1;
            for k = 1:size(bin,1)
                track = tmpdat(tmpdat(:,1) == bin(k),:);
                % percNorm = Track(:,54);
                % percDir = Track(:,55);
                % percConfin = Track(:,56);
                % percAnom = Track(:,57);
                [~, maxIDX] = max(track(:,54:57),[],2);
                if all(maxIDX == maxIDX(1))
                    track(:,end+1) = counter;
                    newTmpDat(tmpdat(:,1) == track(1,1),:) = track;
                    counter = counter+1;
                else
                    %more segments
                    %determine changepoint
                    changepoints = ones(length(maxIDX),1);
                    for m = 1:(length(maxIDX)-1)
                        val1 = maxIDX(m);
                        val2 = maxIDX(m+1);
                        if val1 ~= val2
                            changepoints(m) = 0;
                        end
                    end
                    trackSegsEndIdx = find(changepoints == 0);
                    segNumbers = zeros(size(track,1),1);
                    strt = 1;
                    for j = 1:length(trackSegsEndIdx)
                        segNumbers(strt:trackSegsEndIdx(j)) = counter;
                        strt = trackSegsEndIdx(j)+1;
                        counter = counter+1;
                    end
                    segNumbers(segNumbers == 0) = counter;
                    counter = counter +1;
                    track(:,end+1) = segNumbers;
                    newTmpDat(tmpdat(:,1) == track(1,1),:) = track;
                end
            end
            TrackData{i,1} = newTmpDat;
        end
    else
        %
    end
end