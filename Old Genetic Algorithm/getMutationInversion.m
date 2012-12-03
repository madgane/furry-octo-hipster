function [xChromozome yChromozome] = getMutationInversion(xChromozome,yChromozome,fitnessGain)

probMutation = 1;
if (rand < probMutation)

    betaOne = 1.2;betaTwo = 10;
    pMutation = 1 ./ (betaOne + betaTwo * (std(fitnessGain) / mean(fitnessGain)));

    pMutation = pMutation';
    prbMutationX = rand(size(xChromozome));prbMutationY = rand(size(yChromozome));
    prbMutationX = prbMutationX > pMutation;prbMutationY = prbMutationY > pMutation;

    xChromozome = xor(xChromozome,prbMutationX);
    yChromozome = xor(yChromozome,prbMutationY);
    
end
