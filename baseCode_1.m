
if ~isequal(SimParams.ExtRun,'true')
    clc;clear;
    SimParams.ExtRun = 'false';
end

fwkInitialization;
SimParams.ChannelModel = 'IID';

if ~isequal(SimParams.ExtRun,'true')
    SimParams.SchedType = 'GreedyScheduling';
    SimParams.PrecodingMethod = 'Best_MZF_Method';
end

SimParams.nBands = 1;
SimParams.nBases = 1;
SimParams.userIdices = [10:20:100];

SimParams.PF_dur = 40;
SimParams.maxRank = 4;
SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 1;

if ~isequal(SimParams.ExtRun,'true')
    SimParams.nDrops = 1000;
    SimParams.snrIndex = [25];
end

SimParams.maxUsers = 100;
SimParams.estError = 0.00;
plGlobalOverMaxUsers = -rand(SimParams.nBases,SimParams.maxUsers) * 30;
SimParams.Thrpt = zeros(length(SimParams.snrIndex),length(SimParams.userIdices));

profile on

for iSNR = 1:length(SimParams.snrIndex)
    
    SimParams.iSNR = iSNR;
    SimParams.sPower = 10.^(SimParams.snrIndex(iSNR)/10);
    
    for iUser = 1:length(SimParams.userIdices)
         
        SimParams.nUsers = SimParams.userIdices(1,iUser);
        SimParams.PL_Profile = plGlobalOverMaxUsers(:,1:SimParams.nUsers);
        SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers));
        
        utilityScale = SimParams.nDrops * SimParams.muxRank * SimParams.nBands;
        SimParams.fairness = zeros(length(SimParams.snrIndex),SimParams.nUsers);
        SimStructs.userStruct = cell(SimParams.nUsers,1);SimStructs.baseStruct = cell(SimParams.nBases,1);

        [SimParams,SimStructs] = systemInitialize(SimParams,SimStructs);
        [SimParams,SimStructs] = systemLinking(SimParams,SimStructs);
        
        for iDrop = 1:SimParams.nDrops
            SimParams.iDrop = iDrop;
            [SimParams,SimStructs] = dropInitialize(SimParams,SimStructs);
            [SimParams,SimStructs] = getScheduledUsers(SimParams,SimStructs);
            [SimParams,SimStructs] = getPMatrix(SimParams,SimStructs);
            [SimParams,SimStructs] = performReception(SimParams,SimStructs);
        end

        for cUser = 1:SimParams.nUsers
            SimParams.PFmetric(iSNR,iUser) = SimStructs.userStruct{cUser}.PFmetric;
            SimParams.fairness(iSNR,iUser) = SimStructs.userStruct{cUser}.tAllocation / utilityScale;
            SimParams.Thrpt(iSNR,iUser) = SimParams.Thrpt(iSNR,iUser) + (SimStructs.userStruct{cUser}.crThrpt - 1) / SimParams.nDrops;
        end
    
    end
    
    cState = sprintf('SINR completed - %d',SimParams.snrIndex(iSNR));disp(cState);
    
end

profile off

% profile viewer

SimParams.profiler.schX = SimParams.profiler.schX / (iSNR * iDrop);

if ~isequal(SimParams.ExtRun,'true')
    
    markerS = 'rh-';
    figure(2);hold on;
    plot(SimParams.userIdices,SimParams.Thrpt,markerS);
    xlabel('Number of Users');ylabel('Sum Capacity in Bits/Sec/Hz');grid on;
    hold off;
    
end

% figure(2);hold all;
% plot(SimParams.snrIndex,std(SimParams.Thrpt,0,2));
% xlabel('SNR in dB');ylabel('Capacity Deviation across Users in Bits/Sec/Hz');grid on;

% figure(3);hold all;
% plot(SimParams.snrIndex,std(SimParams.fairness,0,2));
% xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');grid on;
