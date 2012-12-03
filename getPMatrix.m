function [SimParams,SimStructs] = getPMatrix(SimParams,SimStructs)

switch (SimParams.PrecodingMethod)
    
    case 'Best_ZF_Method'
        [SimParams,SimStructs] = getZFMatrix(SimParams,SimStructs);
        
    case 'Best_BF_Method'
        [SimParams,SimStructs] = getBFMatrix(SimParams,SimStructs);
        
    case 'Best_MZF_Method'
        [SimParams,SimStructs] = getMZFMatrix(SimParams,SimStructs);
        
    case 'Best_Network_Method'
        [SimParams,SimStructs] = getNetworkMatrix(SimParams,SimStructs);
        
    case 'Best_DP_Method'
        [SimParams,SimStructs] = getDPMatrix(SimParams,SimStructs);
        
    case 'Best_MaxMinSINR_Method'
        [SimParams,SimStructs] = performMaxMinSINR(SimParams,SimStructs);
        
    case 'Best_WSR_Method'
        [SimParams,SimStructs] = getWSRMatrix(SimParams,SimStructs);
	
    case 'Best_WMMSE_Method'
        [SimParams,SimStructs] = getWMMSEMatrix(SimParams,SimStructs);
        
    case 'Best_TDMZF_Method'
        [SimParams,SimStructs] = getTDMZFMatrix(SimParams,SimStructs);
        
    case 'Best_CZF_Method'
        [SimParams,SimStructs] = getCZFMatrix(SimParams,SimStructs);
        
    otherwise
        display('Unknown Precoder Design !');
        
end
