function [xChromozome yChromozome] = getChildChromosomes(GPStruct)

nParents = 2;
xyGenesP = cell(nParents,1);
avgFitness = GPStruct.xGeneFitnessMean;
prbXselection = GPStruct.xGeneFitness / avgFitness;

rndPRBgen = rand(1,GPStruct.Np);
[~,srtI] = sort(rndPRBgen .* prbXselection,'descend');

for iParent = 1:nParents
    xyGenesP{iParent,1} = GPStruct.geneGroup(:,:,srtI(iParent));
end

xChromozome = xyGenesP{1,1};yChromozome = xyGenesP{2,1};
xChromozome(:,GPStruct.xOverPattern) = xyGenesP{2,1}(:,GPStruct.xOverPattern);
yChromozome(:,GPStruct.xOverPattern) = xyGenesP{1,1}(:,GPStruct.xOverPattern);

end
