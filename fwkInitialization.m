
function [SimParams,SimStructs] = fwkInitialization(SimParams,SimStructs)

stream = RandStream.getGlobalStream;reset(stream);

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.trafficStats.pktService = zeros(length(SimParams.maxArrival),SimParams.nDrops);
end

plModel = char(SimParams.pathLossModel);
uscore_index = find(plModel == '_');

if ~isempty(uscore_index)
    pathLossModel = plModel(1:uscore_index(1,1) - 1);
else
    pathLossModel = plModel;
end    

switch pathLossModel
    
    case 'Isolated'
        
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        for iUser = 1:SimParams.nUsers
            modNode = mod(iUser - 1,SimParams.nBases) + 1;
            SimParams.PL_Profile(:,iUser) = -1e5;
            SimParams.PL_Profile(modNode,iUser) = 0;
        end
        
    case 'CellEdge'
        
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        for iUser = 1:SimParams.nUsers
            modNode = mod(iUser - 1,SimParams.nBases) + 1;
            SimParams.PL_Profile(:,iUser) = -1/1e5;
            SimParams.PL_Profile(modNode,iUser) = 0;
        end
        
    case 'Random'
        
        plGain = str2double(plModel(uscore_index(1,1) + 1:end));
        SimParams.PL_Profile = -rand(SimParams.nBases,SimParams.nUsers) * plGain;
        
    otherwise
        
        xdB = char(SimParams.pathLossModel);xI = strfind(xdB,'_');
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        for iUser = 1:SimParams.nUsers
            modNode = mod(iUser - 1,SimParams.nBases) + 1;
            SimParams.PL_Profile(:,iUser) = str2double(xdB(xI+1:end));
            SimParams.PL_Profile(modNode,iUser) = 0;
        end
end

switch SimParams.arrivalDist    
    case 'Random'        
        pktRatio = 2;
        pktStep = SimParams.maxArrival(1,SimParams.iPkt) / pktRatio;
        randArrival = pktStep:pktStep:pktStep * pktRatio * pktRatio;
    case 'Constant'        
        randArrival = SimParams.maxArrival(1,SimParams.iPkt);        
end

if ~strcmp(SimParams.arrivalDist,'Fixed')
    maxMult = 10;xLength = length(randArrival);
    randIndices = randperm(SimParams.nUsers * maxMult,SimParams.nUsers);
    randIndices = mod(randIndices(1,1:SimParams.nUsers) - 1,xLength) + 1;
else
    randIndices = 1:SimParams.nUsers;
    randArrival = SimParams.FixedPacketArrivals;    
end

SimParams.avgPktValues = randArrival(1,randIndices);

dopplerType = char(SimParams.DopplerType);
uscore_index = find(dopplerType == '_');
maxDoppler = (0.1 * (1 / SimParams.sampTime));

if ~isempty(uscore_index)
    dopType = dopplerType(1:uscore_index(1,1) - 1);
    currDoppler = str2double(dopplerType(uscore_index(1,1) + 1:end));
    if currDoppler > maxDoppler
        currDoppler = maxDoppler;
        display('Warning Doppler Frequency is too high for the current sampling time !');
        display('Resetting the Doppler Frequency !');
    end
else
    dopType = dopplerType;
end    

switch dopType
    case 'Uniform'
        dopplerRealizations = rand(SimParams.nUsers,1) * currDoppler;
    case 'Constant'
        dopplerRealizations = ones(SimParams.nUsers,1) * currDoppler;
end

SimParams.userDoppler = round(dopplerRealizations);
SimStructs.JakesChStruct = cell(SimParams.nUsers,SimParams.nBases,SimParams.nBands);

if strcmp(SimParams.ChannelModel,'Jakes')
    for iUser = 1:SimParams.nUsers
        currentDoppler = SimParams.userDoppler(iUser,1);
        for iBase = 1:SimParams.nBases
            for iBand = 1:SimParams.nBands
                SimStructs.JakesChStruct{iUser,iBase,iBand} = mimochan(SimParams.nTxAntenna,SimParams.nRxAntenna,SimParams.sampTime,currentDoppler);
                SimStructs.JakesChStruct{iUser,iBase,iBand}.ResetBeforeFiltering = 0;
            end
        end
    end
end

nFeedbackOverCycle = 1e6;
SimParams.updateFeedback = zeros(SimParams.nUsers,1);

for iUser = 1:SimParams.nUsers
    SimParams.updateFeedback(iUser,1) = round((1 / SimParams.userDoppler(iUser,1)) / nFeedbackOverCycle / SimParams.sampTime) + 1;
end
