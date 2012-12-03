function [xChromozome yChromozome] = performInversion(xChromozome,yChromozome)

pMutation = 0.5;
[nPhenotypes,geneLength] = size(xChromozome);

if (rand < pMutation)

    aPhenotype = xChromozome;bPhenotype = yChromozome;
    prbMutationX = sort(randi(geneLength,nPhenotypes,2),2);prbMutationY = sort(randi(geneLength,nPhenotypes,2),2);
    
    for iPhenotype = 1:nPhenotypes
        aPhenotype(iPhenotype,prbMutationX(iPhenotype,1):prbMutationX(iPhenotype,2)) = fliplr(xChromozome(iPhenotype,prbMutationX(iPhenotype,1):prbMutationX(iPhenotype,2)));
        bPhenotype(iPhenotype,prbMutationY(iPhenotype,1):prbMutationY(iPhenotype,2)) = fliplr(yChromozome(iPhenotype,prbMutationY(iPhenotype,1):prbMutationY(iPhenotype,2)));
    end   
    
    xChromozome = aPhenotype;yChromozome = bPhenotype;
    
end

end
