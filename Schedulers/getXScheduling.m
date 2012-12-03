function [SimParams,SimStructs] = getXScheduling(SimParams,SimStructs)

schLogic = 'StreamSearch';

for iBand = 1:SimParams.nBands
    
    switch schLogic
        
        case 'EqualShare'
            
            for iBase = 1:SimParams.nBases
                
                cUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
                eH = SimStructs.linkChan{iBase,iBand}(:,:,cUsers);
                augE = reshape(eH(:),SimParams.nTxAntenna,length(cUsers));
                
                [~,~,sortA] = qr(augE,'vector');
                muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
                schedUsers = cUsers(sortA(1,1:muxIFFreeRank),1);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(muxIFFreeRank,1);
                
            end
            
        case 'GroupSS'
            
            nIterations = 10;
            muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
            activeUsers = cell(SimParams.nBases,1);
            
            for iIter = 1:nIterations
                
                subPspace = cell(SimParams.nBases,1);
                
                for iBase = 1:SimParams.nBases
                    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    actSet = uIndices(activeUsers{iBase,1});
                    for jBase = 1:SimParams.nBases
                        if iBase ~= jBase
                            eH = SimStructs.linkChan{jBase,iBand}(:,:,actSet);
                            X = reshape(eH(:),SimParams.nTxAntenna,length(actSet));
                            subPspace{jBase,1} = [subPspace{jBase,1} X];
                        end
                    end
                end
                
                Aspace = cell(SimParams.nBases,1);
                
                for iBase = 1:SimParams.nBases
                    
                    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);kUsers = length(uIndices);
                    X = reshape(eH(:),SimParams.nTxAntenna * SimParams.nRxAntenna,kUsers);
                    
                    
                    ppVolume = zeros(1,kUsers);
                    
                    
                    for iStream = 1:muxIFFreeRank
                        for iUser = 1:kUsers
                            U = [subPspace{iBase,1} , Aspace{iBase,1} , X(:,iUser)];
                            ppVolume(1,iUser) = real(det(U' * U));
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');
                        activeUsers{iBase,1}(1,iStream) = sortI(1,1);
                        Aspace{iBase,1} = X(:,activeUsers{iBase,1}(1,1:iStream));
                    end
                end
                
            end
            
            for iBase = 1:SimParams.nBases
                uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                schedUsers = uIndices(activeUsers{iBase,1},1);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(muxIFFreeRank,1);
            end
            
        case 'InstaSS'
            
            nIterations = 10;
            subPspace = cell(SimParams.nBases,1);
            muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
            activeUsers = cell(SimParams.nBases,1);
            
            for iIter = 1:nIterations
                
                Aspace = cell(SimParams.nBases,1);
                
                for iBase = 1:SimParams.nBases
                    
                    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);kUsers = length(uIndices);
                    X = reshape(eH(:),SimParams.nTxAntenna * SimParams.nRxAntenna,kUsers);
                    
                    ppVolume = zeros(1,kUsers);
                    for iStream = 1:muxIFFreeRank
                        for iUser = 1:kUsers
                            U = [subPspace{iBase,1} , Aspace{iBase,1} , X(:,iUser)];
                            ppVolume(1,iUser) = real(det(U' * U));
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');
                        activeUsers{iBase,1}(1,iStream) = sortI(1,1);
                        Aspace{iBase,1} = X(:,activeUsers{iBase,1}(1,1:iStream));
                    end
                    
                    subPspace = cell(SimParams.nBases,1);
                    for cBase = 1:SimParams.nBases
                        mIndices = SimStructs.baseStruct{cBase,1}.linkedUsers;
                        actSet = mIndices(activeUsers{cBase,1},1);
                        for jBase = 1:SimParams.nBases
                            if cBase ~= jBase
                                eH = SimStructs.linkChan{jBase,iBand}(:,:,actSet);
                                X = reshape(eH(:),SimParams.nTxAntenna,length(actSet));
                                subPspace{jBase,1} = [subPspace{jBase,1}  X];
                            end
                        end
                    end
                end
                
            end
            
            for iBase = 1:SimParams.nBases
                uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                schedUsers = uIndices(activeUsers{iBase,1},1);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(muxIFFreeRank,1);
            end
            
        case 'OrthoSS'
            
            nIterations = 10;
            subPspace = cell(SimParams.nBases,1);
            subDspace = cell(SimParams.nBases,1);
            muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
            activeUsers = cell(SimParams.nBases,1);
            
            for iIter = 1:nIterations
                
                Aspace = cell(SimParams.nBases,1);
                Dspace = cell(SimParams.nBases,1);
                
                for mBase = 1:SimParams.nBases
                    if isempty(subPspace{mBase,1})
                        subPspace{mBase,1} = zeros(SimParams.nTxAntenna);
                    end
                    if isempty(subDspace{mBase,1})
                        subDspace{mBase,1} = zeros(SimParams.nTxAntenna);
                    end
                    if isempty(Aspace{mBase,1})
                        Aspace{mBase,1} = zeros(SimParams.nTxAntenna);
                    end
                    if isempty(Dspace{mBase,1})
                        Dspace{mBase,1} = zeros(SimParams.nTxAntenna);
                    end
                end
                
                
                for iBase = 1:SimParams.nBases
                    
                    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);kUsers = length(uIndices);
                    X = reshape(eH(:),SimParams.nTxAntenna * SimParams.nRxAntenna,kUsers);
                    
                    ppVolumeN = zeros(1,kUsers);
                    ppVolumeD = zeros(1,kUsers);
                    for iStream = 1:muxIFFreeRank
                        for iUser = 1:kUsers
                            U = subPspace{iBase,1} + Aspace{iBase,1};
                            K = subDspace{iBase,1} + Dspace{iBase,1};
                            if ~isempty(U)
                                ppVolumeN(1,iUser) = norm(U' * X(:,iUser));
                            else
                                ppVolumeN(1,iUser) = norm(X(:,iUser));
                            end
                            if ~isempty(K)
                                ppVolumeD(1,iUser) = norm(K' * X(:,iUser));
                            else
                                ppVolumeD(1,iUser) = 1;
                            end
                        end
                        
                        if iStream ~= 1
                            ppVolume = atan(ppVolumeN ./ ppVolumeD);
                        else
                            ppVolume = ppVolumeN;
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');
                        activeUsers{iBase,1}(1,iStream) = sortI(1,1);
                        M = X(:,activeUsers{iBase,1}(1,1:iStream));
                        Aspace{iBase,1} = eye(SimParams.nTxAntenna) - M * inv(M' * M) * M';
                        Dspace{iBase,1} = M * inv(M' * M) * M';
                    end
                    
                    subPspace = cell(SimParams.nBases,1);
                    for mBase = 1:SimParams.nBases
                        if isempty(subPspace{mBase,1})
                            subPspace{mBase,1} = zeros(SimParams.nTxAntenna);
                        end
                        if isempty(subDspace{mBase,1})
                            subDspace{mBase,1} = zeros(SimParams.nTxAntenna);
                        end
                        if isempty(Aspace{mBase,1})
                            Aspace{mBase,1} = zeros(SimParams.nTxAntenna);
                        end
                        if isempty(Dspace{mBase,1})
                            Dspace{mBase,1} = zeros(SimParams.nTxAntenna);
                        end
                    end
                    
                    
                    for cBase = 1:SimParams.nBases
                        mIndices = SimStructs.baseStruct{cBase,1}.linkedUsers;
                        actSet = mIndices(activeUsers{cBase,1},1);
                        for jBase = 1:SimParams.nBases
                            if cBase ~= jBase
                                eH = SimStructs.linkChan{jBase,iBand}(:,:,actSet);
                                X = reshape(eH(:),SimParams.nTxAntenna,length(actSet));
                                subPspace{jBase,1} = subPspace{jBase,1} + (eye(SimParams.nTxAntenna) - X * inv(X' * X) * X');
                                subDspace{jBase,1} = subPspace{jBase,1} + X * inv(X' * X) * X';
                            end
                        end
                    end
                end
                
            end
            
            for iBase = 1:SimParams.nBases
                uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                schedUsers = uIndices(activeUsers{iBase,1},1);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(muxIFFreeRank,1);
            end
            
        case 'StreamSearch'
            
            nIterations = 10;
            muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
            
            Dspace = cell(SimParams.nBases,1);
            Nspace = cell(SimParams.nBases,1);   
            activeUsers = cell(SimParams.nBases,1);
            activeStreams = cell(SimParams.nBases,1);
            completeH = cell(SimParams.nUsers,SimParams.nBases);
            
            for iUser = 1:SimParams.nUsers
                
                bNode = SimStructs.userStruct{iUser,1}.baseNode;
                Hd = SimStructs.linkChan{bNode,iBand}(:,:,iUser);
                [U,~,~] = svd(Hd);
                
                for iBase = 1:SimParams.nBases
                    Hk = SimStructs.linkChan{iBase,iBand}(:,:,iUser);
                    M = U' * Hk;M = M(1:SimParams.maxRank,:);
                    completeH{iUser,iBase} = [completeH{iUser,iBase} M.'];
                end
                
            end
            
            for iIter = 1:nIterations
                
                for iBase = 1:SimParams.nBases
                    
                    Dspace{iBase,1} = [];
                    lkUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    ppVolume = zeros(1,length(lkUsers) * SimParams.maxRank);
                    for iStream = 1:muxIFFreeRank
                        for iUser = 1:length(lkUsers)
                            for jStream = 1:SimParams.maxRank
                                stIndex = (iUser - 1) * SimParams.maxRank + jStream;
                                U = [Nspace{iBase,1} Dspace{iBase,1} completeH{lkUsers(iUser,1),iBase}(:,jStream)];
                                ppVolume(1,stIndex) = real(det(U' * U));
                            end                            
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');
                        activeUsers{iBase,1}(1,iStream) = lkUsers(floor((sortI(1,1) - 1)/SimParams.maxRank) + 1,1);
                        activeStreams{iBase,1}(1,iStream) = mod((sortI(1,1) - 1),SimParams.maxRank) + 1;
                        Dspace{iBase,1} = [Dspace{iBase,1} completeH{activeUsers{iBase,1}(1,iStream),iBase}(:,activeStreams{iBase,1}(1,iStream))];                        
                    end
                    
                    Nspace = cell(SimParams.nBases,1);
                    for kBase = 1:SimParams.nBases
                        if kBase ~= iBase
                            for iUser = 1:length(activeUsers{iBase,1})
                                cUser = activeUsers{iBase,1}(1,iUser);cStream = activeStreams{iBase,1}(1,iUser);
                                Nspace{kBase,1} = [Nspace{kBase,1} completeH{cUser,kBase}(:,cStream)];
                            end
                        end
                    end                    
                end                
            end
            
            for iBase = 1:SimParams.nBases
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = activeUsers{iBase,1}';
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = activeStreams{iBase,1}';
            end
            
    end
    
end
