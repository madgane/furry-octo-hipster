function [SimParams,SimStructs] = performReception(SimParams,SimStructs)

singleNode = 1;
for iBand = 1:SimParams.nBands
    
    linkChannel = SimStructs.actualChannel;
    
    for iUser = 1:SimParams.nUsers
        
        userActive = 0;
        performCooperation = 0;
        cUser = SimStructs.userStruct{iUser,1};
        
        baseNode = cUser.baseNode;
        neighNode = cUser.neighNode;
        
        if length(baseNode) ~= singleNode
            performCooperation = 1;
        end
        
        for kBase = 1:length(baseNode)
            if ~isempty(find(iUser == SimStructs.baseStruct{baseNode(1,kBase),1}.assignedUsers{iBand,1}))
                userActive = 1;break;
            end
        end
        
        if userActive
            I = SimParams.N * eye(SimParams.nRxAntenna);
        else
            I = 1;
        end
        
        S = 0;
        Nacc = 0;
        
        if userActive
            
            for iBase = 1:length(baseNode)
                
                baseIndex = baseNode(1,iBase);
                cBase = SimStructs.baseStruct{baseIndex,1};
                
                gP = cBase.P{iBand,1};
                pIndices = find(iUser == cBase.assignedUsers{iBand,1});
                
                P = gP(:,pIndices);
                H = linkChannel{baseIndex,iBand}(:,:,iUser);
                S = H * P + S;
                
                % Intra Stream Calculation
                
                pIndices = find(iUser ~= cBase.assignedUsers{iBand,1});
                P = gP(:,pIndices);

                if ~isempty(P)
                    N = H * P;
                    if performCooperation
                        Nacc = Nacc + N;
                    else
                        I = I + N * N';
                    end
                end
                
            end
            
            % Inter Stream Calculation
            
            for iBase = 1:length(neighNode)
                
                baseIndex = neighNode(1,iBase);
                cBase = SimStructs.baseStruct{baseIndex,1};
                
                gP = cBase.P{iBand,1};
                pIndices = find(iUser ~= cBase.assignedUsers{iBand,1});
                P = gP(:,pIndices);
                
                H = linkChannel{baseIndex,iBand}(:,:,iUser);
                
                if ~isempty(P)
                    N = H * P;
                    I = I + N * N';
                end
                
            end
            
            if performCooperation
                I = I + Nacc * Nacc';
            end
            
        end
        
        L = eye(size(I)) + (S * S') / (I);
        xThrpt = log2(real(det(L)));
        
        if ~isnan(xThrpt)
            SimStructs.userStruct{iUser,1}.crThrpt = SimStructs.userStruct{iUser,1}.crThrpt + xThrpt;
            SimStructs.userStruct{iUser,1}.lastThrpt = xThrpt + SimStructs.userStruct{iUser,1}.lastThrpt;
        end
        
        if userActive 
            if sign(xThrpt)
                SimStructs.userStruct{iUser,1}.tAllocation = SimStructs.userStruct{iUser,1}.tAllocation + 1;
            end
        end        
        
    end
    
end

% aH = SimStructs.linkChan;
% for iBand = 1:SimParams.nBands
%
%     for iBase = 1:SimParams.nBases
%         pickUsers = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
%
%         txUsers = unique(pickUsers);
%         gP = SimStructs.baseStruct{iBase,1}.P{iBand,1};
%
%         for iUser = 1:length(txUsers)
%
%             cUser = txUsers(iUser,1);
%             stLocs = find(cUser == pickUsers);
%             W = SimStructs.userStruct{cUser,1}.W{iBand,1};
%             H = aH{iBase,iBand}(:,:,cUser);
%
%             W = (W * W');
%             S = H * gP(:,stLocs);
%             I = SimParams.N * eye(size(W));
%
%             for jUser = 1:length(txUsers)
%                 mUser = txUsers(jUser,1);
%                 if mUser ~= cUser
%                     itLocs = find(cUser ~= pickUsers);
%                     N = H * gP(:,itLocs);
%                     I = I + N * N';
%                 end
%             end
%
%             neighNodes = SimStructs.userStruct{cUser,1}.neighNode;
%
%             for jBase = 1:length(neighNodes)
%                 ifJBase = neighNodes(1,jBase);
%                 H = aH{ifJBase,iBand}(:,:,cUser);
%                 P = SimStructs.baseStruct{ifJBase,1}.P{iBand,1};
%                 N = H * P;
%                 I = I + N * N';
%             end
%
%             L = eye(size(I)) + (S * S') / (I);
%             xThrpt = log2(real(det(L)));
%
%             if ~isnan(xThrpt)
%                 SimStructs.userStruct{cUser,1}.lastThrpt = xThrpt;
%                 SimStructs.userStruct{cUser,1}.crThrpt = SimStructs.userStruct{cUser,1}.crThrpt + xThrpt;
%             end
%
%             SimStructs.userStruct{cUser,1}.tAllocation = SimStructs.userStruct{cUser,1}.tAllocation + 1;
%
%         end
%
%     end
%
% end

% for iBand = 1:SimParams.nBands
%
%     for iBase = 1:SimParams.nBases
%
%         pickUsers = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
%         pickStreams = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
%
%         for iRank = 1:length(pickUsers)
%             cStream = pickStreams(iRank,1);cUser = pickUsers(iRank,1);
%             W = SimStructs.userStruct{cUser,1}.W{iBand,1}(:,cStream);
%             P = SimStructs.baseStruct{iBase,1}.P{iBand}(:,iRank);
%             H = SimStructs.linkChanCorrupted{iBase,iBand}(:,:,cUser);
%
%             NPower = 0;
%             SPower = (W' * H * P)' * (W' * H * P);
%             for ifRank = 1:length(pickUsers)
%                 if (iRank ~= ifRank)
%                     P = SimStructs.baseStruct{iBase,1}.P{iBand}(:,ifRank);
%                     NPower = NPower + (W' * H * P)' * (W' * H * P);
%                 end
%             end
%
%             for jBase = 1:SimParams.nBases
%                 if jBase ~= iBase
%                     Hif = SimStructs.linkChanCorrupted{jBase,iBand}(:,:,cUser);
%                     pickUsersJ = SimStructs.baseStruct{jBase,1}.assignedUsers{iBand,1};
%                     for jRank = 1:length(pickUsersJ)
%                         P = SimStructs.baseStruct{jBase,1}.P{iBand}(:,jRank);
%                         NPower = NPower + (W' * Hif * P)' * (W' * Hif * P);
%                     end
%                 end
%             end
%
%             xThrpt = log2(1 + ((SPower) / (trace(W * W') * SimParams.N + NPower)));
%
%             if ~isnan(xThrpt)
%                 SimStructs.userStruct{cUser,1}.lastThrpt = xThrpt;
%                 SimStructs.userStruct{cUser,1}.crThrpt = SimStructs.userStruct{cUser,1}.crThrpt + xThrpt;
%             end
%
%             SimStructs.userStruct{cUser,1}.tAllocation = SimStructs.userStruct{cUser,1}.tAllocation + 1;
%         end
%     end
% end
