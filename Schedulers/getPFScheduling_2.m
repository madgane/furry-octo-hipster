function [SimParams,SimStructs] = getMPFScheduling_2(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    iR = zeros(kUsers,1);
    pF = zeros(kUsers,1);
    mGains = zeros(kUsers,1);
    
    for iUser = 1:kUsers
        pF(iUser,1) = SimStructs.userStruct{uIndices(iUser,1)}.PFmetric;
    end

    augV = zeros(SimParams.nTxAntenna,kUsers);
    for iBand = 1:SimParams.nBands
        
        eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
        for iUser = 1:kUsers
                       
            [~,D,V] = svd(eH(:,:,iUser));
            
            augV(:,iUser) = V(:,1);
            mGains(iUser,1) = 1;
            iR(iUser,1) = log2(1 + D(1,1).^2 * SimParams.sPower) / pF(iUser,1);
            
        end
      
        pairingMatrix = getPairingMatrix(SimParams,iR,augV,mGains);
        xStream = min(SimParams.muxRank,length(pairingMatrix));
        SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = uIndices(pairingMatrix);
        SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(xStream,1);

    end
    
end

end


