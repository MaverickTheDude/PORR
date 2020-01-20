function [dY, dq, dp, sigma] = HDCAfun4R(t, Y, input)
% 4-link pendulum in joint coordinates
p = Y(1:4);     H = [0; 0; 1]; 
q = Y(5:8);     D = [eye(2); [0 0]];
data = getS(q, input);

% assembly
ksi = ksiCoefs(q, p, input);
Q1 = setForcesAtH1(q, input);
assemblyA = struct('Q1',Q1(:,1),'ksi',ksi(1),'S12',data(1).S12);
assemblyB = struct('Q1',Q1(:,2),'ksi',ksi(2),'S12',data(2).S12);
assemblyAB = assemblyPhase(assemblyA, assemblyB);
assemblyA = struct('Q1',Q1(:,3),'ksi',ksi(3),'S12',data(3).S12);
assemblyB = struct('Q1',Q1(:,4),'ksi',ksi(4),'S12',data(4).S12);
assemblyCD = assemblyPhase(assemblyA, assemblyB);
[Q1Sacc, ksiS, S12S] = assemblyPhase(assemblyAB, assemblyCD);

% base body connection
c = - D.' * ksiS.i11 * D; % no inverse here (!)
T1S = D * (c \ D.') * ksiS.i10;     % wzor (43)
T2S = zeros(3,1);
Q1Sart = Q1Sacc; % articulated forces
Q2Sart = zeros(3,1);
V1S = ksiS.i11 * T1S + ksiS.i10;    % wzor (42)
V2S = ksiS.i21 * T1S + ksiS.i20;    % wzor (38)

% disassembly
structS = struct('ksi',ksiS,'T1',T1S,'T2',T2S,'Q1',Q1Sart,'Q2',Q2Sart);
[structCtmp, structDtmp] = disassemblyPhase(structS);

[structA, structB] = disassemblyPhase(structCtmp);
[structC, structD] = disassemblyPhase(structDtmp);

Q1Aart = structA.Q1;   T1A = structA.T1;   T2A = structA.T2;
Q1Bart = structB.Q1;   T1B = structB.T1;   T2B = structB.T2;
Q1Cart = structC.Q1;   T1C = structC.T1;   T2C = structC.T2;
Q1Dart = structD.Q1;   T1D = structD.T1;   T2D = structD.T2;
sigma = [T1A(1:2); T1B(1:2); T1C(1:2); T1D(1:2)];

% Velocity calculation
% wzor (46) zgadza sie ze wzorami (15), (16).
[V1A, V2A] = calculateVelocities(ksi(1), T1A, T2A); % wzory (15), (16)
[V1B, V2B] = calculateVelocities(ksi(2), T1B, T2B);
[V1C, V2C] = calculateVelocities(ksi(3), T1C, T2C);
[V1D, V2D] = calculateVelocities(ksi(4), T1D, T2D);
dq = zeros(4,1);
dq(1) = H' * V1A;
dq(2) = H' * (V1B - V2A);
dq(3) = H' * (V1C - V2B);
dq(4) = H' * (V1D - V2C);

% Articulated momenta
P1Aart = T1A + H*p(1);
P1Bart = T1B + H*p(2);
P1Cart = T1C + H*p(3);
P1Dart = T1D + H*p(4);
P2Aart = T2A - H*p(2);

dS = getdS(q, dq, input);
dSC1_A = dS(1).dSc1;   dSC2_A = dS(1).dSc2;
dSC1_B = dS(2).dSc1;   dSC2_B = dS(2).dSc2;
dSC1_C = dS(3).dSc1;   dSC2_C = dS(3).dSc2;
dSC1_D = dS(4).dSc1;   dSC2_D = dS(4).dSc2;

des3 = H'*(dSC2_C - dSC1_D)*P1Dart;
des2 = des3 + H'*(dSC2_B - dSC1_C)*P1Cart;
des1 = des2 + H'*(dSC2_A - dSC1_B)*P1Bart;

dp = zeros(4,1);
dp(1) = H' * (Q1Aart - dSC1_A * P1Aart) + des1;
dp(2) = H' * (Q1Bart - dSC1_B * P1Bart) + des2;
dp(3) = H' * (Q1Cart - dSC1_C * P1Cart) + des3;
dp(4) = H' * (Q1Dart - dSC1_D * P1Dart);

dY = [dp; dq];
