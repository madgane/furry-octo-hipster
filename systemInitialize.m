

function [SimParams,SimStructs] = systemInitialize(SimParams,SimStructs)

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.crThrpt = 1;
    SimStructs.userStruct{iUser,1}.PFmetric = 0;
    SimStructs.userStruct{iUser,1}.lastThrpt = 0;
    SimStructs.userStruct{iUser,1}.userID = iUser;
    SimStructs.userStruct{iUser,1}.tAllocation = 0;
    SimStructs.userStruct{iUser,1}.W = cell(SimParams.nBands,1);
    
    SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt = 0;
    SimStructs.userStruct{iUser,1}.trafficConfig.bufferLength = 'Inf';
    SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime = zeros(1,SimParams.nDrops);
end

for iBase = 1:SimParams.nBases
    SimStructs.baseStruct{iBase,1}.profile = 0;
    SimStructs.baseStruct{iBase,1}.baseID = iBase;
    SimStructs.baseStruct{iBase,1}.allocGains = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase,1}.allocPattern = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase,1}.P = cell(SimParams.nBands,1);
end

SimParams.profiler.schX = zeros(SimParams.nBases,1);
[SimParams,SimStructs] = generateUserTrafficArrivals(SimParams,SimStructs);

SimStructs.linkChan = cell(SimParams.nBases,SimParams.nBands);
SimStructs.actualChannel = cell(SimParams.nBases,SimParams.nBands);

end