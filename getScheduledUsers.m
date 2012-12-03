function [SimParams,SimStructs] = getScheduledUsers(SimParams,SimStructs)

baseTime = zeros(SimParams.nBases,1);
for iBase = SimParams.nBases
    baseTime(iBase,1) = cputime;
end

chScheduler = char(SimParams.SchedType);
uScoreIndex = find(chScheduler == '_');
if isempty(uScoreIndex)
    scheduleMethod = SimParams.SchedType;
else
    scheduleMethod = chScheduler(1:uScoreIndex(1,1) - 1);
end

switch scheduleMethod
    
    case 'RRScheduling'
        [SimParams,SimStructs] = getRoundRobinScheduling(SimParams,SimStructs);
        
    case 'RouletteScheduling'
        [SimParams,SimStructs] = getRouletteWheelScheduling(SimParams,SimStructs);
    
    case 'GreedyScheduling'
        [SimParams,SimStructs] = getGreedyScheduling(SimParams,SimStructs);
        
    case 'PFScheduling'
        [SimParams,SimStructs] = getFairnessScheduling(SimParams,SimStructs);
        
    case 'BDScheduling'
        [SimParams,SimStructs] = getBDScheduling(SimParams,SimStructs);
        
    case 'PFBDScheduling'
        [SimParams,SimStructs] = getPFBDScheduling(SimParams,SimStructs);
        
    case 'GeneticScheduling'
        [SimParams,SimStructs] = getGeneticScheduling_4(SimParams,SimStructs);
        
    case 'ExScheduling'
        [SimParams,SimStructs] = getExhaustiveScheduling(SimParams,SimStructs);        
        
    case 'XScheduling'
        [SimParams,SimStructs] = getXScheduling(SimParams,SimStructs);   
        
    case 'NetworkScheduling'
        [SimParams,SimStructs] = getNetworkScheduling(SimParams,SimStructs);
        
    case 'CoordScheduling'
        [SimParams,SimStructs] = getCoordinateScheduling(SimParams,SimStructs);
        
    otherwise
        display('Unknown Scheduling Type');
end

for iBase = SimParams.nBases
    SimParams.profiler.schX(iBase,1) = (cputime - baseTime(iBase,1)) + SimParams.profiler.schX(iBase,1);
end

