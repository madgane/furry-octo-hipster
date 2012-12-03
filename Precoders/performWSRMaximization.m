
% nUsers = 6;
% nReceive = 1;
% nTransmit = 4;
% 
% W = cell(nUsers,1);
% H = cell(nUsers,1);
% for iUser = 1:nUsers
%     H{iUser,1} = complex(randn(nReceive,nTransmit),randn(nReceive,nTransmit)) / sqrt(2);    
% end
% 
% WSRMaximization_t.nLayers = nUsers;
% WSRMaximization_t.pFactor = ones(nUsers,1);
% WSRMaximization_t.H = H;
% WSRMaximization_t.Pt = 100;

function [WSRMaximization_t] = performWSRMaximization(WSRMaximization_t)

% Initalizing P and W to SVD variables.

oldNorm = 5e10;
WSRMaximization_t.cvx_status = 'Failed';
PwrVector = ones(WSRMaximization_t.nLayers,1);
[~,nTransmit] = size(WSRMaximization_t.H{1,1});
P = zeros(nTransmit,WSRMaximization_t.nLayers);

for iLayer = 1:WSRMaximization_t.nLayers
    [U D V] = svd(WSRMaximization_t.H{iLayer,1});
    P(:,iLayer) = V(:,1);
end

WSRMaximization_t.P = P;
for iLayer = 1:WSRMaximization_t.nLayers
    WSRMaximization_t.W{iLayer,1} = calculateMMSEWvector(WSRMaximization_t.H,WSRMaximization_t.P,iLayer);
end

SNRconstraints = 1e-5 * ones(WSRMaximization_t.nLayers,1);

while (1)
    
    [WSRMaximization_t,~,SNRconstraints_n] = performSGNProg(WSRMaximization_t,SNRconstraints);    
    [WSRMaximization_t] = performSOCPprog(WSRMaximization_t,SNRconstraints_n);
    
    SNRgamma = sprintf('Gamma - %s',mat2str(single(SNRconstraints_n')));
    WSRstatus = sprintf('SGN / SOCP - %s / %s',WSRMaximization_t.cvx_status.sgn,WSRMaximization_t.cvx_status.socp);
    disp(WSRstatus);disp(SNRgamma);

    if strcmp(WSRMaximization_t.cvx_status.sgn,'Solved')
        if strcmp(WSRMaximization_t.cvx_status.socp,'Solved')
            newNorm = norm((SNRconstraints_n - SNRconstraints) ./ SNRconstraints);
            if (oldNorm - newNorm) < 1e-2
                break;
            else
                oldNorm = newNorm;
            end
        end
    end
    
    SNRconstraints = SNRconstraints_n;
    for iLayer = 1:WSRMaximization_t.nLayers
        PwrVector(iLayer,1) = norm(WSRMaximization_t.P(:,iLayer)).^2;
        WSRMaximization_t.P(:,iLayer) = WSRMaximization_t.P(:,iLayer) / norm(WSRMaximization_t.P(:,iLayer));
    end
    
    PwrVector = PwrVector * (WSRMaximization_t.Pt / sum(PwrVector));
    
    for iLayer = 1:WSRMaximization_t.nLayers
        WSRMaximization_t.W{iLayer,1} = calculateMMSEWvector(WSRMaximization_t.H,WSRMaximization_t.P,iLayer);
        SNRconstraints(iLayer,1) = calculateSINR(WSRMaximization_t.H,WSRMaximization_t.P,WSRMaximization_t.W{iLayer,1},iLayer,PwrVector);
    end
   
end

