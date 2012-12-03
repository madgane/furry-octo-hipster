function [xChromozome yChromozome] = performMutation(xChromozome,yChromozome,GPStruct)

probMutation = 0.5;
if (rand < probMutation)

    betaOne = 1.2;betaTwo = 10;
    pMutation = 1 ./ (betaOne + betaTwo * (GPStruct.xGeneFitnessDeviation / GPStruct.xGeneFitnessMean));

    prbMutationX = rand(size(xChromozome));prbMutationY = rand(size(yChromozome));
    prbMutationX = prbMutationX > pMutation;prbMutationY = prbMutationY > pMutation;

    xChromozome = xor(xChromozome,prbMutationX);
    yChromozome = xor(yChromozome,prbMutationY);
    
end

end
