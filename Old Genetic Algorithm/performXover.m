function [xChromozome yChromozome] = performXover(fitnessValues,xGenes)

nParents = 2;
[nUsers,nPop] = size(xGenes);

genePicks = 1:2:nUsers;
xyGenesP = cell(nParents,1);
avgFitness = mean(fitnessValues);
prbXselection = fitnessValues / avgFitness;

rndPRBgen = rand(1,nPop);
[~,srtI] = sort(rndPRBgen .* prbXselection,'descend');

for iParent = 1:nParents
    xyGenesP{iParent,1} = xGenes(:,srtI(iParent));
end

xChromozome = xyGenesP{1,1};yChromozome = xyGenesP{2,1};
xChromozome(genePicks,1) = xyGenesP{2,1}(genePicks);
yChromozome(genePicks,1) = xyGenesP{1,1}(genePicks,1);

end
