function [dSc1, dSc2] = getdS2(fi, dfi, input)
Om = [0, -1; 1, 0];
M = input.M;                sC2_loc = input.sC2_loc;
sC1_loc = input.sC1_loc;    %s12_loc = sC1_loc-sC2_loc;
[dSc1, dSc2] = deal(zeros(3));

sC1 = Rot(fi)*sC1_loc;
sC2 = Rot(fi)*sC2_loc;

dSc1(3,1:2) = (-dfi*sC1).'; 
dSc2(3,1:2) = (-dfi*sC2).';