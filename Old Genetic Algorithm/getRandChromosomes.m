function [XY] = getRandChromosomes(E,K,Nt,Np)

XY = rand(E,K,Np);
XY = 1 .* (XY >= (K - Nt)/K);

for iNp = 1:Np
    for iE = 1:E
        if (sum(XY(iE,:,iNp)) > Nt)
            while (1)
                xOne = find(XY(iE,:,iNp) == 1);
                pPickIndex = floor(rand(1,1) * length(xOne)) + 1;
                XY(iE,xOne(pPickIndex),iNp) = 0;
                if (sum(XY(iE,:,iNp)) == Nt)
                    break;
                end
            end
        end
        
        if (sum(XY(iE,:,iNp)) == 0)
            while (1)
                xOne = find(XY(iE,:,iNp) == 0);
                pPickIndex = floor(rand(1,1) * length(xOne)) + 1;
                XY(iE,xOne(pPickIndex),iNp) = 1;
                if (sum(XY(iE,:,iNp)) ~= 0)
                    break;
                end
            end
        end
    end
end

end

  