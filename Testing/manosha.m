
clc;
clear all;

K = 100;
rand('seed',0);randn('seed',0);

PL = rand(100,1) * -50 - 120;

NV = -114;

TxPwr = 43;


SINR = TxPwr + PL - NV

NV = 0;

TxPwr = 43 + 114;


SINR = TxPwr + PL - NV