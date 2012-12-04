
function [SimParams SimStructs] = getCNetworkBFWMMSEDesign(SimParams,SimStructs)

iIter = 0;
maxIter = 1e4;
epsilonCheck = min(1e-4,(SimParams.sPower)^(-2));
nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);

perBSpowerConstraint = 'false';

if strcmp(perBSpowerConstraint,'true')
    
    for iBand = 1:SimParams.nBands
        
        continueAgain = 1;
        W = cell(SimParams.nUsers,1);
        U = cell(SimParams.nUsers,1);
        V = cell(SimParams.nUsers,1);
        Haug = cell(SimParams.nUsers,1);
        linkChannel = SimStructs.linkChan;
        
        for iUser = 1:SimParams.nUsers
            for iBase = 1:SimParams.nBases
                Haug{iUser,1} = [Haug{iUser,1} linkChannel{iBase,iBand}(:,:,iUser)];
            end
            
            V{iUser,1} = complex(ones(SimParams.nTxAntenna * SimParams.nBases,1:SimParams.maxRank),ones(SimParams.nTxAntenna * SimParams.nBases,1:SimParams.maxRank));
            V{iUser,1} = sqrt(SimParams.sPower / (SimParams.nUsers / SimParams.nBases)) * V{iUser,1} / trace(V{iUser,1}' * V{iUser,1});
        end
        
        while continueAgain
            
            % U Matrix calculation
            
            for iUser = 1:SimParams.nUsers
                
                H = Haug{iUser,1};
                J = eye(SimParams.nRxAntenna) * SimParams.N;
                
                for jUser = 1:SimParams.nUsers
                    HV = H * V{jUser,1};
                    J = J + HV * HV';
                end
                
                U{iUser,1} = inv(J) * H * V{iUser,1};
                W{iUser,1} = inv(eye(nStreams) - U{iUser,1}' * H * V{iUser,1});
            end
            
            for iBase = 1:SimParams.nBases
                
                Isum = 0;Dsum = 0;
                for iUser = 1:SimParams.nUsers
                    cUser = SimStructs.userStruct{iUser,1};
                    H_HU = linkChannel{iBase,iBand}(:,:,iUser)' * U{iUser,1};
                    Isum = Isum + cUser.weighingFactor * H_HU * W{iUser,1} * H_HU';
                    W_2 = W{iUser,1} * W{iUser,1};
                    Dsum = Dsum + cUser.weighingFactor^2 * H_HU * W_2 * H_HU';
                end
                
                mu_star = bisectionEstimateMU(Isum,Dsum,SimParams.sPower);
                
                Isum = Isum + mu_star * eye(SimParams.nTxAntenna);
                sI = (iBase - 1) * SimParams.nTxAntenna + 1;eI = sI + SimParams.nTxAntenna - 1;
                
                Iinv = pinv(Isum);
                for iUser = 1:SimParams.nUsers
                    cUser = SimStructs.userStruct{iUser,1};
                    Vb = cUser.weighingFactor * Iinv * linkChannel{iBase,iBand}(:,:,iUser)' * U{iUser,1} * W{iUser,1};
                    V{iUser,1}(sI:eI,:) = Vb;
                end
                
            end
            
            clc;
            reshape(cell2mat(V),SimParams.nBases * SimParams.nTxAntenna,SimParams.nUsers)
            
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
        end
        
        % Assigning the V and U to the corresponding users
        
        linkUsers = 1:SimParams.nUsers;
        for iBase = 1:SimParams.nBases
            
            P = [];aUsers = [];
            sI = (iBase - 1) * SimParams.nTxAntenna + 1;
            eI = sI + SimParams.nTxAntenna - 1;
            SimStructs.baseStruct{iBase,1}.linkedUsers = linkUsers';
            
            for iUser = 1:SimParams.nUsers
                P = [P V{iUser,1}(sI:eI,:)];
                aUsers = [aUsers ; repmat(iUser,nStreams,1)];
            end
            
            SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = aUsers;
            
        end
        
        for iUser = 1:SimParams.nUsers
            SimStructs.userStruct{iUser,1}.W{iBand,1} = U{iUser,1};
            SimStructs.userStruct{iUser,1}.baseNode = 1:SimParams.nBases;
            SimStructs.userStruct{iUser,1}.neighNode = [];
        end
        
    end
    
else
    
    for iBand = 1:SimParams.nBands
        
        continueAgain = 1;
        W = cell(SimParams.nUsers,1);
        U = cell(SimParams.nUsers,1);
        V = cell(SimParams.nUsers,1);
        Haug = cell(SimParams.nUsers,1);
        linkChannel = SimStructs.linkChan;
        
        for iUser = 1:SimParams.nUsers
            for iBase = 1:SimParams.nBases
                Haug{iUser,1} = [Haug{iUser,1} linkChannel{iBase,iBand}(:,:,iUser)];
            end
            
            V{iUser,1} = complex(ones(SimParams.nTxAntenna * SimParams.nBases,1:SimParams.maxRank),ones(SimParams.nTxAntenna * SimParams.nBases,1:SimParams.maxRank));
            V{iUser,1} = sqrt(SimParams.sPower / (SimParams.nUsers / SimParams.nBases)) * V{iUser,1} / trace(V{iUser,1}' * V{iUser,1});
        end
        
        while continueAgain
            
            % U Matrix calculation
            
            for iUser = 1:SimParams.nUsers
                
                H = Haug{iUser,1};
                J = eye(SimParams.nRxAntenna) * SimParams.N;
                
                for jUser = 1:SimParams.nUsers
                    HV = H * V{jUser,1};
                    J = J + HV * HV';
                end
                
                U{iUser,1} = inv(J) * H * V{iUser,1};
                W{iUser,1} = inv(eye(nStreams) - U{iUser,1}' * H * V{iUser,1});
            end
            
            Isum = 0;Dsum = 0;
            for iUser = 1:SimParams.nUsers
                cUser = SimStructs.userStruct{iUser,1};
                H_HU = Haug{iUser,1}' * U{iUser,1};
                Isum = Isum + cUser.weighingFactor * H_HU * W{iUser,1} * H_HU';
                W_2 = W{iUser,1} * W{iUser,1};
                Dsum = Dsum + cUser.weighingFactor * H_HU * W_2 * H_HU';
            end
            
            mu_star = bisectionEstimateMU(Isum,Dsum,SimParams.sPower * SimParams.nBases);
            Isum = Isum + mu_star * eye(SimParams.nTxAntenna * SimParams.nBases);
            
            Iinv = pinv(Isum);
            for iUser = 1:SimParams.nUsers
                cUser = SimStructs.userStruct{iUser,1};
                V{iUser,1} = cUser.weighingFactor * Iinv * Haug{iUser,1}' * U{iUser,1} * W{iUser,1};
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
        end
        
        % Assigning the V and U to the corresponding users
        
        linkUsers = 1:SimParams.nUsers;
        for iBase = 1:SimParams.nBases
            
            P = [];aUsers = [];
            sI = (iBase - 1) * SimParams.nTxAntenna + 1;
            eI = sI + SimParams.nTxAntenna - 1;
            SimStructs.baseStruct{iBase,1}.linkedUsers = linkUsers';
            
            for iUser = 1:SimParams.nUsers
                P = [P V{iUser,1}(sI:eI,:)];
                aUsers = [aUsers ; repmat(iUser,nStreams,1)];
            end
            
            SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = aUsers;
            
        end
        
        for iUser = 1:SimParams.nUsers
            SimStructs.userStruct{iUser,1}.W{iBand,1} = U{iUser,1};
            SimStructs.userStruct{iUser,1}.baseNode = 1:SimParams.nBases;
            SimStructs.userStruct{iUser,1}.neighNode = [];
        end
        
    end
    
end
