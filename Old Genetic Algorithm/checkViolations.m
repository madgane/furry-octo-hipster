function [xChild,yChild] = checkViolations(xChild,yChild,SimParams)

maxN = SimParams.muxRank;nUsers = length(xChild);

while 1
    if isempty(find(sum(xChild) == (1:maxN)))
        if sum(xChild) == 0
            xChild(randi(nUsers,1,1),1) = 1;
        else
            xLocs = find(xChild == 1);
            xChild(xLocs(randi(length(xLocs),1,1)),1) = 0;
        end
    else
        break;
    end
end

while 1
    if isempty(find(sum(yChild) == (1:maxN)))
        if sum(yChild) == 0
            yChild(randi(nUsers,1,1),1) = 1;
        else
            yLocs = find(yChild == 1);
            yChild(yLocs(randi(length(yLocs),1,1)),1) = 0;
        end
    else
        break;
    end
end

