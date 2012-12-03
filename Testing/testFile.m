

clc;

Bo = randi([4,10],1,100000);
Wo = randi([1,8],1,100000);

B = Bo;W = Wo;

P = mean(abs(B - W));
D = (W + P) > B;
OverallcostMean = sum(~D) * P

P = median(abs(B - W));
D = (W + P) > B;
OverallcostMedian = sum(~D) * P

runCost = [];

runInt = 1;
for i = 0:1:100
    D = (W + i) > B;
    runCost(runInt,1) = sum(~D) * i;
    runInt = runInt + 1;
end

hold all;
plot(runCost)
plot(repmat(OverallcostMean,1,length(runCost)));
plot(repmat(OverallcostMedian,1,length(runCost)));
