function xGeneFitness = getChromosomeWeights(SimParams,SimStructs,GPStruct)

E = SimParams.nBases;
xGeneFitness = zeros(1,GPStruct.Np);

for iPool = 1:GPStruct.Np
    
    cUserIndices = cell(E,1);
    cGene = GPStruct.geneGroup(:,:,iPool);
    
    for iBase = 1:E
        cUserIndices{iBase,1} = find(cGene(iBase,:) == 1);
    end
    
    switch SimParams.PrecodingMethod
        
        case 'Best_ZF_Method'
            xWeight = getCompleteZFprecoders(SimParams,SimStructs,cUserIndices,GPStruct.cBand);
            
        case 'Best_WMMSE_Method'
            
    end
    
    for iBase = 1:E
        for iUser = 1:length(cUserIndices{iBase,1})
            xGeneFitness(1,iPool) = xGeneFitness(1,iPool) + ...
                xWeight{iBase,1}(1,iUser) / SimStructs.userStruct{cUserIndices{iBase,1}(1,iUser)}.PFmetric;
        end            
    end
    
end
