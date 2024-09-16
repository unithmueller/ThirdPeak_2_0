flatData = [];
                    dataNames = [];
                    gloIterator = 1;
                    for i = 1:size(gaussDat{1,1})
                        for j = 1:size(gaussDat{1,1}{i,1})
                            dataNames(gloIterator) = gaussDat{1,2}{i,2}(j,1);
                            flatData(gloIterator,1:2) = [gaussDat{1,1}{i,1}(j),gaussDat{1,1}{i,2}(j)];
                            gloIterator = gloIterator+1;
                        end
                    end