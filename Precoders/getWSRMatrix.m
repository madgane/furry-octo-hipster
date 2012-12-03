function [SimParams,SimStructs] = getWSRMatrix(SimParams,SimStructs)

WSRMaximization_t.Pt = SimParams.sPower;
WSRMaximization_t.nLayers = SimParams.muxRank;
WSRMaximization_t.H = cell(SimParams.muxRank,1);
WSRMaximization_t.W = cell(SimParams.muxRank,1);
WSRMaximization_t.pFactor = ones(SimParams.muxRank,1);

for iBase = 1:SimParams.nBases
    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    for iBand = 1:SimParams.nBands
        
        W = cell(kUsers,1);
        eG = zeros(SimParams.maxRank,kUsers);
        eH = SimStructs.linkChan(:,:,uIndices,iBase,iBand);
        for iUser = 1:kUsers
            d = svd(eH(:,:,iUser));
            eG(1:length(d),iUser) = d;
            [W{iUser,1}, ~, ~] = svd(eH(:,:,iUser));
        end
        
        [~, sortIndex] = sort(eG(:),'descend');sortIndex = sortIndex - 1;
        
        for iRank = 1:SimParams.maxRank
            if (iRank > SimParams.muxRank)
                break;
            end
            cUser = floor(sortIndex(iRank,1) / SimParams.maxRank) + 1;
            cStream = mod(sortIndex(iRank,1),SimParams.maxRank) + 1;
            
            WSRMaximization_t.H{iRank,1} = eH(:,:,cUser);
            WSRMaximization_t.W{iRank,1} = W{cUser,1}(:,cStream);
        end
        
        [WSRMaximization_t] = performWSRMaximization(WSRMaximization_t);
        SimStructs.baseStruct{iBase}.P(:,:,iBand) = WSRMaximization_t.P;
        SimStructs.baseStruct{iBase}.allocPattern(:,iBand) = sortIndex;
        SimStructs.baseStruct{iBase}.allocGains(:,:,iBand) = eG;
        
        for iRank = 1:SimParams.muxRank
            cUser = floor(sortIndex(iRank,1) / SimParams.maxRank) + 1;
            cStream = mod(sortIndex(iRank,1),SimParams.maxRank) + 1;            
            SimStructs.userStruct{uIndices(cUser,1)}.W(:,cStream,iBand) = WSRMaximization_t.W{iRank,1};
        end
        
    end
    
end

end
