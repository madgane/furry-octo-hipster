
clc;clear;
pathAddition;
SimParams.DebugMode = 'false';

SimParams.ExtRun = 'false';
SimParams.allActive = 'false';

SimParams.ChannelModel = 'IID';
SimParams.pathLossModel = 'Random_10';
SimParams.DopplerType = 'Uniform_100';

SimParams.weighingEqual = 'false';
SimParams.SchedType = 'BDScheduling_SP';
SimParams.PrecodingMethod = 'Best_ZF_Method';
SimParams.weightedSumRateMethod = 'StreamScheduling';

SimParams.nDrops = 50;
SimParams.snrIndex = [15];

SimParams.PF_dur = 40;
SimParams.sampTime = 1e-3;
SimParams.estError = 0.00;

SimParams.nBands = 1;
SimParams.nBases = 1;
SimParams.nUsers = 50;

SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 1;

SimParams.gracePeriod = 0;
SimParams.arrivalDist = 'Constant';

SimParams.maxArrival = 5;
SimParams.FixedPacketArrivals = [50,50,50,50,50,50,5,5,5,5];

nSINRSamples = length(SimParams.snrIndex);
nPacketSamples = length(SimParams.maxArrival);
SimStructs.userStruct = cell(SimParams.nUsers,1);
SimStructs.baseStruct = cell(SimParams.nBases,1);

SimParams.maxRank = min(SimParams.nRxAntenna,SimParams.nTxAntenna);
SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers / SimParams.nBases));

SimParams.Thrpt = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
utilityScale = SimParams.nDrops * SimParams.muxRank * SimParams.nBands;
SimParams.fairness = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);

queueBacklogs = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
queueBacklogsOverTime = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples,SimParams.nDrops);

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

        cState = sprintf('SINR completed - %d',SimParams.snrIndex(iSNR));disp(cState);
    end
    
end

testPktIndex = 1;testSINRIndex = 1;
SimParams.profiler.schX = SimParams.profiler.schX / (iSNR * iDrop);
   
fairnessS = squeeze(SimParams.fairness(:,:,testPktIndex));
throughputSumS = squeeze(SimParams.Thrpt(:,:,testPktIndex));

queueBackLogsS = squeeze(queueBacklogs(testSINRIndex,:,:));
queueBacklogsOverTimeS = squeeze(queueBacklogsOverTime(testSINRIndex,:,testPktIndex,:));

if strcmp(SimParams.allActive,'true')
    
    markerS = 'o-';
    figure(1);hold all;
    SimParams.sumThrpt = sum(throughputSumS,2);
    plot(SimParams.snrIndex,SimParams.sumThrpt,markerS);
    xlabel('SNR in dB');ylabel('Sum Capacity in Bits/Sec/Hz');grid on;

%     figure(2);hold all;
%     JainMean = mean(throughputSumS,2).^2;JainVar = var(throughputSumS,0,2);
%     JainIndex_capacity = JainMean ./ (JainMean + JainVar);
%     
%     plot(SimParams.snrIndex,JainIndex_capacity,markerS);
%     xlabel('SNR in dB');ylabel('Capacity Deviation across Users in Bits/Sec/Hz');grid on;
    
%     figure(3);hold all;
%     JainMean = mean(fairnessS,2).^2;JainVar = var(fairnessS,0,2);
%     JainIndex_utility = JainMean ./ (JainMean + JainVar);
%     
%     plot(SimParams.snrIndex,JainIndex_utility,markerS);
%     xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');grid on;

else
    
    markerS = '.-';
    figure(5);hold all;
    plot(1:SimParams.nDrops,sum(queueBacklogsOverTimeS,1));
    xlabel('Slot Index');ylabel('Queue Backlogs (pkts) over Time');grid on;
    
end
