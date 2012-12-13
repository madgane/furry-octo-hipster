
function [SimParams,SimStructs] = fwkInitialization(SimParams,SimStructs)

stream = RandStream.getGlobalStream;reset(stream);

% Path Loss Model related code

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
        
    case 'Fixed'
        
        SimParams.PL_Profile = SimParams.PL_Profile;
        
    otherwise
        
        xdB = char(SimParams.pathLossModel);xI = strfind(xdB,'_');
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        for iUser = 1:SimParams.nUsers
            modNode = mod(iUser - 1,SimParams.nBases) + 1;
            SimParams.PL_Profile(:,iUser) = str2double(xdB(xI+1:end));
            SimParams.PL_Profile(modNode,iUser) = 0;
        end
end

% Queue Related code

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.trafficStats.pktService = zeros(length(SimParams.maxArrival),SimParams.nDrops);
end

switch SimParams.arrivalDist    
    case 'Random'        
        nRandomness = 100;
        randArrival = rand(1,nRandomness) * SimParams.maxArrival(1,SimParams.iPkt);
    case 'Constant'        
        randArrival = SimParams.maxArrival(1,SimParams.iPkt);        
end

if ~strcmp(SimParams.arrivalDist,'Fixed')
    randIndices = randi([1,length(randArrival)],1,SimParams.nUsers);
else
    randIndices = 1:SimParams.nUsers;
    randArrival = SimParams.FixedPacketArrivals;
end

SimParams.avgPktValues = randArrival(1,randIndices);
[SimParams,SimStructs] = generateUserTrafficArrivals(SimParams,SimStructs);

% Doppler / Small scale related code

legacychannelsim(true);

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

SimParams.updateFeedback = zeros(SimParams.nUsers,1);

for iUser = 1:SimParams.nUsers
    feedbackCycle = (1 / SimParams.userDoppler(iUser,1)) * SimParams.fbFraction;
    SimParams.updateFeedback(iUser,1) = round(feedbackCycle / SimParams.sampTime) + 1;
end

end
