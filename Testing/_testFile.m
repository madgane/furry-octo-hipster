

clc;

B = [2 3 3 2];
W = [1 2 2 1];

P = mean(B) - mean(W);

D = (W + P) > B;

Overallcost = W .* ~D + sum(~D) * P + B .* D;

