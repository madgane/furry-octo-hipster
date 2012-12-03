function pairingMatrix = getPairingMatrix(SimParams,iR,augV,mGains)

Cx = abs(augV' * augV);
Cx = Cx * diag(mGains);

Ng = 10;Np = 10;
vecGenes = randi([0 1],SimParams.nUsers,Np);
vecGenes = checkConstraintViolations(vecGenes,SimParams);

newGenes = zeros(size(vecGenes));
fitnessGain = calcFitness(vecGenes,Cx);


for iG = 1:Ng
    
    for iChild = 1:2:(Np - 2)
        
        [xChild,yChild] = performXover(fitnessGain,vecGenes);
        [xChild,yChild] = getMutationInversion(xChild,yChild,fitnessGain);
        
        [xChild,yChild] = checkViolations(xChild,yChild,SimParams);
        newGenes(:,iChild) = xChild;newGenes(:,iChild + 1) = yChild;
        
    end
    
    [~,srtI] = sort(fitnessGain,'descend');
    newGenes(:,iChild + 2) = vecGenes(:,srtI(1,1));
    newGenes(:,iChild + 3) = vecGenes(:,srtI(1,2));
    
    vecGenes = newGenes;
    fitnessGain = calcFitness(vecGenes,Cx);
        
end
    
[~,maxI] = max(fitnessGain);
pairingMatrix = find(vecGenes(:,maxI) == 1);
    
    
end
