function xWeight = getCompleteZFprecoders(SimParams,SimStructs,cUserIndices,iBand)

Hx = SimStructs.linkChan;
P = cell(SimParams.nBases,1);
xWeight = cell(SimParams.nBases,1);
W = cell(SimParams.nBases,SimParams.nTxAntenna);

for iBase = 1:SimParams.nBases
    H = [];
    for iUser = 1:length(cUserIndices{iBase,1})
        [U, ~, ~] = svd(Hx{iBase,iBand}(:,:,cUserIndices{iBase,1}(1,iUser)));
        
        W{iBase,iUser} = U(:,1);
        H = [H ; W{iBase,iUser}(:,1)' * Hx{iBase,iBand}(:,:,cUserIndices{iBase,1}(1,iUser))];
    end
    
    X = pinv(H);
    P{iBase,1} = performWFAlgorithm(X,SimParams.sPower);
end

for iBase = 1:SimParams.nBases
    for iUser = 1:length(cUserIndices{iBase,1})
        Spower = norm(W{iBase,iUser}' * Hx{iBase,iBand}(:,:,cUserIndices{iBase,1}(1,iUser)) * P{iBase,1}(:,iUser))^2;
        Npower = 0;
        for kUser = 1:length(cUserIndices{iBase,1})
            if (iUser ~= kUser)
                Npower = Npower + norm(W{iBase,iUser}' * Hx{iBase,iBand}(:,:,cUserIndices{iBase,1}(1,iUser)) * P{iBase,1}(:,kUser))^2;
            end
        end
        
        for jBase = 1:SimParams.nBases
            if jBase ~= iBase
                for lUser = 1:length(cUserIndices{jBase,1})
                    Npower = Npower + norm(W{iBase,iUser}' * Hx{jBase,iBand}(:,:,cUserIndices{iBase,1}(1,iUser)) * P{jBase,1}(:,lUser))^2;
                end
            end
        end       
        
        xGain = Spower / (Npower + norm(W{iBase,iUser})^2);
        
        xThrpt = log2(1 + xGain);
        xWeight{iBase,1} = [xWeight{iBase,1} , xThrpt];
        
    end
end


