
function [Pmatrix Pavg] = performWFAlgoCVX(PinvMatrix,Pt,Gains)

[Nt,Nu] = size(PinvMatrix);
Pg = diag(PinvMatrix' * PinvMatrix);

cvx_begin gp
    
    variable Pwr(Nu,1)

    maximize sum(log(1 + Pwr./Pg))
    
    subject to
    
    sum(Pwr) <= Pt;
        
cvx_end

Pavg = Pwr;
Pmatrix = PinvMatrix * diag(sqrt(Pwr));


