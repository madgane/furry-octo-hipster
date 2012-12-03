function vecGenes = checkConstraintViolations(vecGenes,SimParams)

[Nusers Npop] = size(vecGenes);

xSum = sum(vecGenes);
for iSum = 1:Npop
    while 1
        if ~length(find(xSum(1,iSum) == [1:SimParams.muxRank]))
            vecGenes(:,iSum) = randi([0 1],Nusers,1);
            xSum(1,iSum) = sum(vecGenes(:,iSum));
        else
            break;
        end
    end        
end
