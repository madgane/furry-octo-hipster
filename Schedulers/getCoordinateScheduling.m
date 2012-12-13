
function [SimParams,SimStructs] = getCoordinateScheduling(SimParams,SimStructs)

assIndex = cell(SimParams.nBases,SimParams.nBands);
assUsers = cell(SimParams.nBases,SimParams.nBands);
assStreams = cell(SimParams.nBases,SimParams.nBands);

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.baseNode = [];
    SimStructs.userStruct{iUser,1}.neighNode = [];
end

for iBand = 1:SimParams.nBands
    
    Haug = [];iIndex = 0;
    uLocs = zeros(SimParams.nUsers * SimParams.maxRank,3);
    
    for iUser = 1:SimParams.nUsers
        for iBase = 1:SimParams.nBases
            H = SimStructs.linkChan{iBase,iBand}(:,:,iUser);
            [U,~,~] = svd(H);M = U' * H;
            if SimParams.queueWt
                M = M.' * (SimStructs.userStruct{iUser,1}.weighingFactor);
            else
                M = M.' * sign(SimStructs.userStruct{iUser,1}.weighingFactor);
            end
            Haug = [Haug M];
            
            for iStream = 1:SimParams.maxRank
                iIndex = iIndex + 1;
                uLocs(iIndex,:) = [iUser iStream iBase];
            end
            
        end
    end
    
    Hm = Haug;
    ppVolume = zeros(iIndex,1);
    for iStream = 1:iIndex
        ppVolume(iStream,1) = real(det(Hm(:,iStream)' * Hm(:,iStream)));
    end
    
    ppVolume = zeros(iIndex,1);
    eMatrix = cell(SimParams.nBases,1);    
    for iRank = 1:SimParams.muxRank * SimParams.nBases;
        for iUser = 1:iIndex            
            cBase = uLocs(iUser,3);
            if isempty(assUsers{cBase,iBand})
                NullSpace = eye(SimParams.nTxAntenna);
            else
                NullSpace = eMatrix{cBase,1};
            end
            M = Hm(:,iUser)' * NullSpace;
            ppVolume(iUser,1) = real(trace(M' * M));
        end
        
        [~,maxI] = max(ppVolume);
        cUser = uLocs(maxI,1);cStream = uLocs(maxI,2);cBase = uLocs(maxI,3);
        
        assIndex{cBase,iBand} = [assIndex{cBase,iBand} ; maxI];
        assUsers{cBase,iBand} = [assUsers{cBase,iBand} ; cUser];
        assStreams{cBase,iBand} = [assStreams{cBase,iBand} ; cStream];        
        eMatrix{cBase,1} = null(Haug(:,assIndex{cBase,iBand}) * Haug(:,assIndex{cBase,iBand})');
        
        userLocs = find(cUser == uLocs(:,1));
        modLocs = uLocs(userLocs,:);
        oLocs = userLocs(cBase ~= modLocs(:,3));
        Hm(:,oLocs) = zeros(SimParams.nTxAntenna,length(oLocs));
                    
    end
    
end

for iBase = 1:SimParams.nBases
   
    usersLinked = [];
    for iBand = 1:SimParams.nBands
        usersLinked = [usersLinked ; assUsers{iBase,iBand}];        
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = assUsers{iBase,iBand};
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = assStreams{iBase,iBand};
    end
    
    SimStructs.baseStruct{iBase,1}.linkedUsers = unique(usersLinked);
end

iBand = 1;
ovBaseNodes = 1:SimParams.nBases;
for iUser = 1:SimParams.nUsers
    assNode = [];
    for iBase = 1:SimParams.nBases
        if ~isempty(find(iUser == assUsers{iBase,iBand}))
            assNode = [assNode iBase];
        end
    end

    if ~isempty(assNode)
        SimStructs.userStruct{iUser,1}.baseNode = assNode;
        SimStructs.userStruct{iUser,1}.neighNode = ovBaseNodes(1,(assNode ~= ovBaseNodes));
    else
        SimStructs.userStruct{iUser,1}.baseNode = [];
        SimStructs.userStruct{iUser,1}.neighNode = [];
    end
end
