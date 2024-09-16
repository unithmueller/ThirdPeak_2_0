function [flatData, dataNames] = unpackDataCell(dataCell)
flatData = [];
dataNames = strings;
gloIterator = 1;
for i = 1:size(dataCell{1,1},1)
    for j = 1:size(dataCell{1,1}{i,1},1)
        dataNames(gloIterator) = dataCell{1,2}{i,2}(j,1);
        flatData(gloIterator,1:2) = [dataCell{1,1}{i,1}(j),dataCell{1,1}{i,2}(j)];
        gloIterator = gloIterator+1;
    end
end
end