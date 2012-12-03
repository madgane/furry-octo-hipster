function [SimParams,SimStructs] = getPFBDScheduling(SimParams,SimStructs)

charScheduling = char(SimParams.SchedType);
uscore_index = find(charScheduling == '_');
caseStudy = charScheduling(uscore_index + 1:end);

for iBand = 1:SimParams.nBands
    
    for iBase = 1:SimParams.nBases
    
        schedUsers = zeros(SimParams.muxRank,1);
        schedStreams = zeros(SimParams.muxRank,1);        
        uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
        Haug = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
        
        iIndex = 0;
        kUsers = length(uIndices);X = [];
        xLocs = zeros(SimParams.maxRank * kUsers,2);
        
        for iUser = 1:kUsers            
            cUser = uIndices(iUser,1);
            [U,~,~] = svd(Haug(:,:,iUser));
            M = U' * Haug(:,:,iUser) * sign(SimStructs.userStruct{cUser,1}.weighingFactor);
            
            for iRank = 1:SimParams.maxRank
                iIndex = iIndex + 1;
                X = [X M(iRank,:).'];
                xLocs(iIndex,:) = [cUser,iRank];                
            end
        end        
        
        switch caseStudy
            
            case 'AHP'
                
                augM = [];
                for iStream = 1:SimParams.muxRank
                    for iUser = 1:kUsers * SimParams.maxRank
                        for jUser = 1:kUsers * SimParams.maxRank
                            if iUser ~= jUser
                                U = [augM X(:,iUser) X(:,jUser)];
                                G(jUser,iUser) = real(det(U' * U));
                            else
                                U = [augM X(:,iUser)];
                                G(iUser,iUser) = real(det(U' * U));
                            end
                        end
                    end
                    
                    [~,~,V] = svd(G);
                    [~,iSort] = sort(abs(V(:,1)),'descend');
                    
                    augM = [augM X(:,iSort(1,1))];
                    schedUsers(iStream,1) = xLocs(iSort(1,1),1);
                    schedStreams(iStream,1) = xLocs(iSort(1,1),2);
                    
                end
                
        end
                
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = schedUsers;
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
        
    end
    
end
