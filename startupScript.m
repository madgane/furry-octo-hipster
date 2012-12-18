
clc;clear;
pathAddition;

SimParams.DebugMode = 'false';
SimParams.queueMode = 'false';

SimParams.ChannelModel = 'IID';
SimParams.pathLossModel = 'Random_20';
SimParams.DopplerType = 'Constant_100';

SimParams.queueWt = 0;
SimParams.weighingEqual = 'true';
SimParams.SchedType = 'XScheduling_InstaSS';
SimParams.PrecodingMethod = 'Best_CZF_Method';
SimParams.weightedSumRateMethod = 'StreamScheduling';

SimParams.nDrops = 10;
SimParams.snrIndex = [-5:5:15];

SimParams.PF_dur = 40;
SimParams.sampTime = 1e-3;
SimParams.estError = 0.00;
SimParams.fbFraction = 0.0;

SimParams.nBands = 1;
SimParams.nBases = 2;
SimParams.nUsers = 20;

SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 1;

SimParams.gracePeriod = 0;
SimParams.arrivalDist = 'Constant';

SimParams.maxArrival = [100];
SimParams.FixedPacketArrivals = [10,10,10,10,10,10,1,1,1,1];
SimParams.PL_Profile = [5 -inf 5 -inf 5 -inf 1e-20 0; -inf 5 -inf 5 -inf 5 0 1e-20];

nSINRSamples = length(SimParams.snrIndex);
nPacketSamples = length(SimParams.maxArrival);
SimStructs.userStruct = cell(SimParams.nUsers,1);
SimStructs.baseStruct = cell(SimParams.nBases,1);

SimParams.maxRank = min(SimParams.nRxAntenna,SimParams.nTxAntenna);
SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers));

SimParams.Thrpt = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
utilityScale = SimParams.nDrops * SimParams.muxRank * SimParams.nBands;
SimParams.fairness = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);

queueBacklogs = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
queueBacklogsOverTime = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples,SimParams.nDrops);
SimParams.txPower = zeros(length(SimParams.maxArrival),length(SimParams.snrIndex),SimParams.nBases);

for iPkt = 1:length(SimParams.maxArrival)
    
    SimParams.iPkt = iPkt;
    [SimParams,SimStructs] = fwkInitialization(SimParams,SimStructs);
    
    for iSNR = 1:length(SimParams.snrIndex)

        SimParams.N = 1;
        SimParams.iSNR = iSNR;
        SimParams.sPower = 10.^(SimParams.snrIndex(iSNR)/10);
        [SimParams,SimStructs] = systemInitialize(SimParams,SimStructs);
        [SimParams,SimStructs] = systemLinking(SimParams,SimStructs);

        % Resetting for every SNRs
        stream = RandStream.getGlobalStream;reset(stream);

        for iDrop = 1:SimParams.nDrops
            SimParams.iDrop = iDrop;
            [SimParams,SimStructs] = dropInitialize(SimParams,SimStructs);
            [SimParams,SimStructs] = getScheduledUsers(SimParams,SimStructs);
            [SimParams,SimStructs] = getPMatrix(SimParams,SimStructs);
            [SimParams,SimStructs] = performReception(SimParams,SimStructs);
        end

        for iUser = 1:SimParams.nUsers
            SimParams.PFmetric(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser}.PFmetric;
            SimParams.fairness(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser}.tAllocation / utilityScale;
            SimParams.Thrpt(iSNR,iUser,iPkt) = (SimStructs.userStruct{iUser}.crThrpt - 1) / (SimParams.nDrops * SimParams.nBands);
            queueBacklogs(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt;
            queueBacklogsOverTime(iSNR,iUser,iPkt,:) = SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime;
        end

        if strcmp(SimParams.DebugMode,'true')
            display(squeeze(queueBacklogs(iSNR,:,iPkt)));
        end
        
        cState = sprintf('SINR completed - %d',SimParams.snrIndex(iSNR));disp(cState);
    end
    
end

SimResults.avgTxPower = SimParams.txPower / SimParams.nDrops;

if strcmp(SimParams.queueMode,'false')
    
    SimResults.sumThrpt = sum(SimParams.Thrpt(:,:,end),2);
    SimResults.thrptFairness = sum(SimParams.fairness(:,:,end),2);
    
    markerS = 'o-';
    figure(1);hold all;
    SimParams.sumThrpt = SimResults.sumThrpt;
    plot(SimParams.snrIndex,SimParams.sumThrpt,markerS);
    xlabel('SNR in dB');ylabel('Sum Capacity in Bits/Sec/Hz');grid on;

%     figure(2);hold all;
%     JainMean = mean(SimResults.sumThrpt,2).^2;JainVar = var(SimResults.sumThrpt,0,2);
%     JainIndex_capacity = JainMean ./ (JainMean + JainVar);
%     
%     plot(SimParams.snrIndex,JainIndex_capacity,markerS);
%     xlabel('SNR in dB');ylabel('Capacity Deviation across Users in Bits/Sec/Hz');grid on;
    
%     figure(3);hold all;
%     JainMean = mean(SimResults.thrptFairness,2).^2;JainVar = var(SimResults.thrptFairness,0,2);
%     JainIndex_utility = JainMean ./ (JainMean + JainVar);
%     
%     plot(SimParams.snrIndex,JainIndex_utility,markerS);
%     xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');grid on;

else
    
    SimResults.queueBackLogs = queueBacklogs;
    SimResults.queueBackLogsOverTime = queueBacklogsOverTime;
    
%     markerS = '.-';
%     figure(5);hold all;
%     plot(1:SimParams.nDrops,sum(squeeze(SimResults.queueBackLogsOverTime(end,:,end,:)),1));
%     xlabel('Slot Index');ylabel('Queue Backlogs (pkts) over Time');grid on;
    
    markerS = '.*-';
    plot(SimParams.maxArrival,sum(squeeze(SimResults.queueBackLogs(end,:,:)),1));
    xlabel('Average Arrival Rate');ylabel('Average Queue Size (pkts)');grid on;
    hold all;
    
end
