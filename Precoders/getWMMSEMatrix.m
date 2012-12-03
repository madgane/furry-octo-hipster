function [SimParams,SimStructs] = getWMMSEMatrix(SimParams,SimStructs)

switch SimParams.weightedSumRateMethod
    
    case 'PreScheduling'
        
        [SimParams,SimStructs] = getPreWeightedMMSEDesign(SimParams,SimStructs);

    case 'PerformScheduling'
        
        [SimParams,SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
        
    case 'DistScheduling'
        
        [SimParams,SimStructs] = getDWeightedMMSEDesign(SimParams,SimStructs);
        
    case 'CNetworkBF'
        
        [SimParams,SimStructs] = getCNetworkBFWMMSEDesign(SimParams,SimStructs);
        
    case 'DNetworkBF'
        
        [SimParams,SimStructs] = getDNetworkBFWMMSEDesign(SimParams,SimStructs);
        
        
    case 'StreamScheduling'
        
        [SimParams,SimStructs] = getStrWeightedMMSEDesign(SimParams,SimStructs);
                
    otherwise
        
        display('Unknown Weighted Sum Rate Option !');
                
end

end
