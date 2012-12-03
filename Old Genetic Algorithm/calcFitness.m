function fitnessGain = calcFitness(vecGenes,Cx)

[~,Npop] = size(vecGenes);
fitnessGain = zeros(1,Npop);

for iPop = 1:Npop
    fitnessM = Cx(vecGenes(:,iPop) == 1,vecGenes(:,iPop) == 1);
    totalFactor = 1 + sum(fitnessM)' - diag(fitnessM);
    tFactor = diag(fitnessM) ./ totalFactor;
    fitnessGain(1,iPop) = sum(log2(1 + tFactor));
end

if length(find(fitnessGain < 0))
    display('Negative Numbers');
    pause;
end

end