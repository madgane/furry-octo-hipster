
function [Pmatrix Pavg] = performQueuedWF(PinvMatrix,Pt,Q)

cvx_quiet('true');
pAllocationType = 'WSR_ZF';

Q = Q / log2(exp(1));
[~,nUsers] = size(PinvMatrix);
pwrUsers = diag(PinvMatrix' * PinvMatrix);

cvx_begin

variable cvxP(nUsers,1)

switch pAllocationType
    
    case 'WSR_ZF'
        
        maximize sum(log(1 + cvxP ./ pwrUsers) .* Q)
        
        subject to
            sum(cvxP) <= Pt;

    case 'DifferentialQueues'
        
        variable t(nUsers,1)        
        minimize(norm(t,2))
        
        subject to
            sum(cvxP) <= Pt;
            max(Q - log(1 + cvxP ./ pwrUsers),0) <= t;
            
end
   
cvx_end

Pavg = cvxP;
Pmatrix = PinvMatrix;
avgPWR = 1./sqrt(pwrUsers);

for iCol = 1:length(avgPWR)
    if (avgPWR(iCol,1) ~= Inf)
        Pmatrix(:,iCol) = Pmatrix(:,iCol) * avgPWR(iCol,1);
    else
        Pmatrix(:,iCol) = Pmatrix(:,iCol) * 0;
    end
end

Pmatrix = Pmatrix * diag(sqrt(Pavg));


