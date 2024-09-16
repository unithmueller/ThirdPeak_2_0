function [finalIDs, trackfilterid] = filterTrackIDsByDataFilterStruc(DataFilterStruct, SaveStructure, concatData, propertyIDs, propertyNames) 
%Function that will determine the remaining track ids based on the filter
%settings set in the DataFilterWindow

    %% grab the filter
    nodes = DataFilterStruct.nodes;
    selectedIds = DataFilterStruct.filterIDS;
    %choose if whole track or track segments should be filtered
    if selectedIds == "Track"
        trackfilterid = 1;
    else
        trackfilterid = 3;
    end
    involvedIDs = cell(size(nodes,1),1);
    %% for every node build a filter
    for i = 1:size(nodes,1)
        %{filtertype, filterprop, filterextraprop, filterminval, filtermaxval, filterlogical}
        data = nodes(i).NodeData;
        propfieldname = data{2};
        extrapropfieldname = split(data{3},".");
        %numerical data to filter for
        if propfieldname ~= "DeepSPT" && propfieldname ~= "BaseProperties"
            if propfieldname == "Classification from External"
                propfieldname = "SwiftParams";
                extrapropfieldname = "type";
            elseif propfieldname == "Diffusion Parameters from External"
                propfieldname = "SwiftParams";
            end
            if size(extrapropfieldname,1) == 1
                fielddata = getfield(SaveStructure,propfieldname, extrapropfieldname{1});
            else
                fielddata = getfield(SaveStructure,propfieldname, extrapropfieldname{1},extrapropfieldname{2});
            end
            %some extra code to manage wierd data
            if propfieldname == "MeanJumpDist"
                tmpid = cell2mat(fielddata{1,trackfilterid});
                tmpdat = fielddata{1,2};
                tmparray = cell(length(tmpid),2);
                for j = 1:size(tmpid)
                    tmparray{j,1} = tmpid(j);
                    tmparray{j,2} = tmpdat(j);
                end
                fielddata = tmparray;
                clear tmpid tmpdat tmparray
            end
            %need to differ between id-value and id-multiValue pairs
            multval = (cell2mat(cellfun(@(x) size(x,1),fielddata,'UniformOutput',false)));
            multval = mean(multval(:,2));
            if multval > 1 %multivalues to ids
                %need to look at every data point
                keptIDs = ones(size(fielddata,1),1).*-1;
                for k = 1:size(fielddata,1)
                    tmp = fielddata(k,:);
                    tmpid = tmp{trackfilterid};
                    tmpdat = tmp{2};
                    if size(tmpdat,2) > 1
                        tmpdat = tmpdat(:,2);
                    end
                    %check for include exclude
                    if data{1} == 1 %include
                        decisionMatrix = tmpdat >= data{4} & tmpdat <= data{5};
                    else
                        decisionMatrix = tmpdat < data{4} | tmpdat > data{5};
                    end
                    % if all points are inculeded, take the track, if one is
                    % wrong, leave it
                    if decisionMatrix
                        keptIDs(k) = tmpid;
                    end
                end
                keptIDs(keptIDs < 0) = [];
                a = 1;
                involvedIDs{i} = {data{6}, keptIDs.'};
                clear fielddata
            else %only one value per id
                %this contains the data we want to filter
                if size(fielddata(1,:),2) == 3
                    dataSizes = cellfun(@(x) size(x,1), fielddata);
                    dataSizes = (dataSizes(:,3));
                    newFieldData = zeros( max(sum(cellfun(@(x) size(x,1), fielddata))) ,3);
                    counter = 1;
                    for k = 1:size(dataSizes,1)
                        newFieldData(counter:counter+dataSizes(k)-1,1) = fielddata{k,1};
                        newFieldData(counter:counter+dataSizes(k)-1,2) = fielddata{k,2};
                        newFieldData(counter:counter+dataSizes(k)-1,3) = fielddata{k,3};
                        counter = counter+dataSizes(k);
                    end
                    
                    fielddata = newFieldData;
                else
                    fielddata = cell2mat(fielddata);
                end
                %check for include exclude
                if data{1} == 1 %include
                    ids = fielddata(fielddata(:,2) >= data{4} & fielddata(:,2) <= data{5},trackfilterid);
                else
                    ids = fielddata(fielddata(:,2) < data{4} | fielddata(:,2) > data{5},trackfilterid);
                end
                ids = unique(ids);
                involvedIDs{i} = {data{6}, ids};
            end
        else
            %it is deepSPT
            if selectedIds ~= "Track"
                trackfilterid = 66;
            end
            propIDIdx = find(propertyNames == extrapropfieldname);
            %check for include exclude
            if data{1} == 1 %include
                idFilter = (concatData(:,propIDIdx) >= data{4}) & (concatData(:,propIDIdx) <= data{5});
            else
                idFilter = (concatData(:,propIDIdx) < data{4} | concatData(:,propIDIdx) > data{5});
            end
            notIdFilter = ~idFilter;
            inIds = unique(concatData(idFilter,trackfilterid));
            outIds = unique(concatData(notIdFilter,trackfilterid));
            mask = ismember(inIds, outIds);
            inIds = inIds(~mask);
            ids = inIds;
            a=1;
            involvedIDs{i} = {data{6}, ids};
        end
    end
    %%  connect the filters by their logicals
    %check that there are more then 1 connection to do
    if size(nodes,1) > 1
        startids = involvedIDs{1}{2};
        for i = 2:size(nodes,1)
            nextids = involvedIDs{i}{2};
            if involvedIDs{i}{1} == 1 %and
                Lia = ismember(startids,nextids);
                startids = startids(Lia);
            else %or
                tmp = [startids; nextids];
                startids = unique(tmp);
            end
        end
        finalIDs = startids;
    else %just one
        finalIDs = involvedIDs{1}{2};
    end
    %% return the final IDs
    finalIDs = sort(finalIDs);
end