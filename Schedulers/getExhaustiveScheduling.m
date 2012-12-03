function [SimParams,SimStructs] = getExhaustiveScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases
    
    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    userArray = getExhaustiveArray(kUsers,SimParams.muxRank);
    [nGroups,~] = size(userArray);
    
    for iBand = 1:SimParams.nBands                
        
        chINVpwr = zeros(nGroups,1);
        eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
        
        for iGroup = 1:nGroups
            
            augH = [];
            for iStream = 1:SimParams.muxRank
                augH = [augH ; eH(:,:,userArray(iGroup,iStream))];
            end
            eP = pinv(augH);
            chINVpwr(iGroup,1) = trace(eP' * eP);
            
        end
        
        [~,minI] = min(chINVpwr);
        SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = uIndices(userArray(minI,:));
        SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(SimParams.muxRank,1);

    end
    
end

end

