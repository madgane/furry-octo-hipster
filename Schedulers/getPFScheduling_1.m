function [SimParams,SimStructs] = getPFScheduling_1(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    pF = zeros(kUsers,1);
    for iUser = 1:kUsers
        pF(iUser,1) = SimStructs.userStruct{uIndices(iUser,1)}.PFmetric;
    end
    
    if ~isequal(SimParams.ExtRun,'true')
        caseStudy = 2;
    else
        caseStudy = SimParams.caseStudy;
    end
    
    switch (caseStudy)
        
        case 1
            
            for iBand = 1:SimParams.nBands
        
                eG = zeros(SimParams.maxRank,kUsers);
                eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
                for iUser = 1:kUsers
                    d = svd(eH(:,:,iUser));
                    eG(1:length(d),iUser) = log2(1 + d.^2 * SimParams.sPower) / pF(iUser,1);
                end

                [~, sortIndex] = sort(eG(:),'descend');sortIndex = sortIndex - 1;
                pickUsers = floor(sortIndex / SimParams.maxRank) + 1;
                pickStreams = mod(sortIndex,SimParams.maxRank) + 1;        
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = uIndices(pickUsers(1:SimParams.muxRank));
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = pickStreams(1:SimParams.muxRank);

            end
            
        case 2
            
            for iBand = 1:SimParams.nBands
        
                X = [];
                pfMetrics = [];
                eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
                
                for iUser = 1:kUsers

                    [U, D, V] = svd(eH(:,:,iUser));
                    X = [X , V * D'];d = diag(D);                   
                    pfMetrics = [pfMetrics ; log2(1 + d(1:SimParams.nRxAntenna).^2 * SimParams.sPower) / pF(iUser,1)];

                end
                
                G = [];
                kUsers = length(pfMetrics);
                quantilePF = median(pfMetrics);
                schedUser = zeros(1,SimParams.muxRank);
                schedLayer = zeros(1,SimParams.muxRank);
                perStreamMetric = zeros(kUsers,1);
                
                for iStream = 1:SimParams.muxRank
                    
                    for iUser = 1:kUsers
                        if pfMetrics(iUser,1) > quantilePF
                            M = [G X(:,iUser)];
                            perStreamMetric(iUser,1) = real(det(M' * M));
                        else
                            perStreamMetric(iUser,1) = 0;
                        end                        
                    end
                    
                    [~,maxI] = max(perStreamMetric);
                    schedUser(1,iStream) = floor(maxI / SimParams.nRxAntenna);
                    schedLayer(1,iStream) = mod(maxI - 1,SimParams.nRxAntenna) + 1;
                    G = [G X(:,maxI)];
                    
                end
                
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = uIndices(schedUser);
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = schedLayer';
                
            end

    end
    
end

end


