
function [Pmatrix Pavg] = performQueuedWF(PinvMatrix,Pt,Q)

[Nt,Nu] = size(PinvMatrix);
Pg = diag(PinvMatrix' * PinvMatrix);

 cvx_quiet('true');

% cvx_begin
%     
%     variable Pwr(Nu,1)
% 
%     maximize sum(log(1 + Pwr./Pg).*Q)
%     
%     subject to
%     
% 
%     sum(Pwr) <= Pt;
%         
% cvx_end

Q = Q / log2(exp(1));

cvx_begin
    
    variable t(Nu,1)
    variable Pwr(Nu,1)

    minimize sum(t)
    
    subject to    

    max(Q - log(1 + Pwr./Pg),0) <= t
    sum(Pwr) <= Pt;
        
cvx_end


Pavg = Pwr;

Pmatrix = PinvMatrix;
avgPWR = (1./sqrt(diag(Pmatrix' * Pmatrix)));
for iCol = 1:length(avgPWR)
    if (avgPWR(iCol,1) ~= Inf)
        Pmatrix(:,iCol) = Pmatrix(:,iCol) * avgPWR(iCol,1);
    else
        Pmatrix(:,iCol) = Pmatrix(:,iCol) * 0;
    end
end

Pmatrix = Pmatrix * diag(sqrt(Pwr));


