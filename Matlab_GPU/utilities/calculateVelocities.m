function [V1, V2] = calculateVelocities(ksi, T1, T2)

V1 = ksi.i11*T1 + ksi.i12*T2 + ksi.i10;
V2 = ksi.i21*T1 + ksi.i22*T2 + ksi.i20;