clc
clear all;

K = 2;
M = 2;

H = [2 -0.5 ; -1 2];
w1 = 1;w2 = 5;

p = zeros(2,1);

while (1)
    
    % step 1
    
    for k = 1:K
        for j = k:K
            sG = 0;
            for l = 1:j
                if l ~= k
                    sG = sG + H(:,l) * H(:,l)' * p(l,1);
                end
            end
            sG = eye(M) + sG;
            alpha(k,j) = H(:,k)' * inv(sG) * H(:,k);
        end
    end
    
    while (1)
        
        pmax = 1e7;pmin = 0;
        pavg = (pmax + pmin) * 0.5;
        
        
        
        
    end
    
    
    p = pn * (1/K) + (1 - (1/K)) * p;
    
    
end