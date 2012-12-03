function xyChromozomes = performConstraintQualification(xyChromozomes,GPStruct)

nXY = length(xyChromozomes);
[E,~] = size(xyChromozomes{1,1});

for iXY = 1:nXY
    for iBase = 1:E
        if sum(xyChromozomes{iXY,1}(iBase,:)) > GPStruct.maxUsers
            while 1
                xPresent = find(xyChromozomes{iXY,1}(iBase,:) == 1);
                xLoc = floor(rand(1,1) * length(xPresent)) + 1;
                xyChromozomes{iXY,1}(iBase,xPresent(xLoc)) = 0;
                
                if (sum(xyChromozomes{iXY,1}(iBase,:)) == GPStruct.maxUsers)
                    break;
                end
            end            
        end
        
        if sum(xyChromozomes{iXY,1}(iBase,:)) == 0
            while 1
                xPresent = find(xyChromozomes{iXY,1}(iBase,:) == 0);
                xLoc = floor(rand(1,1) * length(xPresent)) + 1;
                xyChromozomes{iXY,1}(iBase,xPresent(xLoc)) = 1;
                
                if (sum(xyChromozomes{iXY,1}(iBase,:)) ~= 0)
                    break;
                end
            end            
        end
    end
end

   
