
function [SimParams SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs)

iIter = 0;
maxIter = 1e4;
epsilonCheck = min(1e-4,(SimParams.sPower)^(-2));
nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);

SumCapacity = cell(SimParams.nBands,1);

for iBand = 1:SimParams.nBands
    
    continueAgain = 1;
    W = cell(SimParams.nUsers,1);
    U = cell(SimParams.nUsers,1);
    V = cell(SimParams.nUsers,1);
    linkChannel = SimStructs.linkChan;
    
    Ui = cell(SimParams.nUsers,1);
    
    for iUser = 1:SimParams.nUsers
        V{iUser,1} = complex(ones(SimParams.nTxAntenna,nStreams),ones(SimParams.nTxAntenna,nStreams));
        V{iUser,1} = sqrt(SimParams.sPower / (SimParams.nUsers / SimParams.nBases)) * V{iUser,1} / trace(V{iUser,1}' * V{iUser,1});
    end
    
    while continueAgain
        
        % U Matrix calculation
        
        for iUser = 1:SimParams.nUsers
            
            cUser = SimStructs.userStruct{iUser,1};
            J = eye(SimParams.nRxAntenna) * SimParams.N;
            
            for jUser = 1:SimParams.nUsers
                ifUser = SimStructs.userStruct{jUser,1};
                HV = linkChannel{ifUser.baseNode,iBand}(:,:,iUser) * V{jUser,1};
                J = J + HV * HV';
            end
            
            H = linkChannel{cUser.baseNode,iBand}(:,:,iUser);
            
            HdVd = H * V{iUser,1};
            U{iUser,1} = inv(J) * HdVd;
            W{iUser,1} = inv(eye(nStreams) - U{iUser,1}' * H * V{iUser,1});
            Ui{iUser,1} = J - HdVd * HdVd' - eye(SimParams.nRxAntenna) * SimParams.N;
            
        end
        
        for iBase = 1:SimParams.nBases
            
            cBase = SimStructs.baseStruct{iBase,1};
            linkedUsers = cBase.linkedUsers;
            
            Isum = 0;Dsum = 0;
            for iUser = 1:SimParams.nUsers
                cUser = SimStructs.userStruct{iUser,1};
                H_HU = linkChannel{iBase,iBand}(:,:,iUser)' * U{iUser,1};
                Isum = Isum + cUser.weighingFactor * H_HU * W{iUser,1} * H_HU';
                if cUser.baseNode == iBase
                    W_2 = W{iUser,1} * W{iUser,1};
                    Dsum = Dsum + cUser.weighingFactor^2 * H_HU * W_2 * H_HU';
                end
            end
            
            mu_star = bisectionEstimateMU(Isum,Dsum,SimParams.sPower);
            Isum = Isum + mu_star * eye(SimParams.nTxAntenna);
            
            Iinv = pinv(Isum);
            for iUser = 1:length(linkedUsers)
                cIndex = linkedUsers(iUser,1);
                cUser = SimStructs.userStruct{cIndex,1};
                V{cIndex,1} = cUser.weighingFactor * Iinv * linkChannel{iBase,iBand}(:,:,cIndex)' * U{cIndex,1} * W{cIndex,1};
            end
            
        end
        
        if ~iIter
            continueAgain = 1;
        else
            currDeviation = 0;
            for iUser = 1:SimParams.nUsers
                currDeviation = currDeviation + abs(log(det(W{iUser,1})) - log(det(W_prev{iUser,1})));
            end
            
            if currDeviation < epsilonCheck
                continueAgain = 0;
            end
            
            if iIter > maxIter
                continueAgain = 0;
                display('Lack of Convergence !');
            end
        end
        
        W_prev = W;
        iIter = iIter + 1;
        SumCapacity{iBand,1} = [SumCapacity{iBand,1} ; performMockReception(SimParams,SimStructs,V,iBand)];
        
    end
    
    % Assigning the V and U to the corresponding users
    
    for iBase = 1:SimParams.nBases
        cBase = SimStructs.baseStruct{iBase,1};
        
        cBase.P{iBand,1} = [];
        cBase.assignedUsers{iBand,1} = [];
        cBase.assignedStreams{iBand,1} = [];
        
        for iUser = 1:length(cBase.linkedUsers)
            cUserIndex = cBase.linkedUsers(iUser,1);
            cUser = SimStructs.userStruct{cUserIndex,1};
            cBase.P{iBand,1} = [cBase.P{iBand,1} , V{cUserIndex,1}];
            cUser.W{iBand,1} = U{cUserIndex,1};
            
            xStreams = (1:nStreams)';
            cBase.assignedUsers{iBand,1} = [cBase.assignedUsers{iBand,1} ; repmat(cUserIndex,length(xStreams),1)];
            cBase.assignedStreams{iBand,1} = [cBase.assignedStreams{iBand,1} ; xStreams];
            SimStructs.userStruct{cUserIndex,1} = cUser;
        end
        
        SimStructs.baseStruct{iBase,1} = cBase;
    end
    
end

