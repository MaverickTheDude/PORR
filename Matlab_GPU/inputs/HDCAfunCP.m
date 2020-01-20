function [dY, dq, dp, sigma] = HDCAfunCP(t, Y, input)
% Cart - Pendulum (1 link), joint coordinates
n = length(Y)/2;
p = Y(1:n);   
q = Y(n+1:end);
data = getS(q, input);

% assembly
ksi = ksiCoefs(q, p, input);
ind = findT(t, input);
Q1 = setForcesAtH1(q, input, ind);  % <= HDCA issue
% Gdy zadana sila tlumienia, Q1 trzeba wyznaczyc PO fazie assembly (TO DO)
assemblyA = struct('Q1',Q1(:,1),'ksi',ksi(1),'S12',data(1).S12);
assemblyB = struct('Q1',Q1(:,2),'ksi',ksi(2),'S12',data(2).S12);
[Q1Sacc, ksiS, S12S] = assemblyPhase(assemblyA, assemblyB);
%{
assemblyAB = assemblyPhase(assemblyA, assemblyB);
assemblyC = struct('Q1',Q1(:,3),'ksi',ksi(3),'S12',data(3).S12);
assemblyD = struct('Q1',Q1(:,4),'ksi',ksi(4),'S12',data(4).S12);
assemblyCD = assemblyPhase(assemblyC, assemblyD);
[Q1Sacc, ksiS, S12S] = assemblyPhase(assemblyAB, assemblyCD);
%}

% base body connection
D = input{1}(1).D;
c = - D.' * ksiS.i11 * D; % no inverse here (!)
T1S = D * (c \ D.') * ksiS.i10;     % wzor (43)
T2S = zeros(3,1);
Q1Sart = Q1Sacc; % articulated forces
Q2Sart = zeros(3,1);

% disassembly
structS = struct('ksi',ksiS,'T1',T1S,'T2',T2S,'Q1',Q1Sart,'Q2',Q2Sart);
[structA, structB] = disassemblyPhase(structS);
%{
[structCtmp, structDtmp] = disassemblyPhase(structS);
[structA, structB] = disassemblyPhase(structCtmp);
[structC, structD] = disassemblyPhase(structDtmp);
%}

Q1Aart = structA.Q1;   T1A = structA.T1;   T2A = structA.T2;
Q1Bart = structB.Q1;   T1B = structB.T1;   T2B = structB.T2;
sigma = [T1A(2:3); T1B(1:2)];
%{
Q1Cart = structC.Q1;   T1C = structC.T1;   T2C = structC.T2;
Q1Dart = structD.Q1;   T1D = structD.T1;   T2D = structD.T2;
sigma = [T1A(2:3); T1B(1:2); T1C(1:2); T1D(1:2)];
%}

% Velocity calculation
% wzor (46) zgadza sie ze wzorami (15), (16).
[V1A, V2A] = calculateVelocities(ksi(1), T1A, T2A); % wzory (15), (16)
[V1B, V2B] = calculateVelocities(ksi(2), T1B, T2B);

dq = zeros(2,1);
dq(1) = data(1).H' * V1A;
dq(2) = data(2).H' * (V1B - V2A);
%{
[V1C, V2C] = calculateVelocities(ksi(3), T1C, T2C);
[V1D, V2D] = calculateVelocities(ksi(4), T1D, T2D);
dq(3) = data(3).H' * (V1C - V2B);
dq(4) = data(4).H' * (V1D - V2C);
%}

% Articulated momenta
P1Aart = T1A + data(1).H*p(1);
P1Bart = T1B + data(2).H*p(2);
% P1Cart = T1C + data(3).H*p(3);
% P1Dart = T1D + data(4).H*p(4);
% P2Aart = - P1Bart;
% P2Aart = T2A - data(1).H*p(2);

dS = getdS(q, dq, input);
dSC1_A = dS(1).dSc1;   dSC2_A = dS(1).dSc2;
dSC1_B = dS(2).dSc1;   dSC2_B = dS(2).dSc2;
% dSC1_C = dS(3).dSc1;   dSC2_C = dS(3).dSc2;
% dSC1_D = dS(4).dSc1;   dSC2_D = dS(4).dSc2;

% descendant(i) = H(i+1) * (...)
dp = zeros(2,1);
% {
dp(1) = data(1).H' * (Q1Aart - dSC1_A * P1Aart);
dp(2) = data(2).H' * (Q1Bart - dSC1_B * P1Bart);
%}
%{
des3 = data(4).H'*(dSC2_C - dSC1_D)*P1Dart;
des2 = data(3).H'*(dSC2_B - dSC1_C)*P1Cart + des3;
des1 = 0; % des1 == 0, poniewaz odnosi do p(1), czyli zmiennej zwiazanej z para postepowa
dp(1) = data(1).H' * (Q1Aart - dSC1_A * P1Aart) + des1;
dp(2) = data(2).H' * (Q1Bart - dSC1_B * P1Bart) + des2;
dp(3) = data(3).H' * (Q1Cart - dSC1_C * P1Cart) + des3;
dp(4) = data(4).H' * (Q1Dart - dSC1_D * P1Dart);
%}
dY = [dp; dq];

% Spelnienie rownan (12,13)
%{
M1V1_LHS = data(1).M1 *  V1A; % LHS
M1V1_RHS = P1Aart+data(1).S12*P2Aart;  % RHS [ok]
M2V2_LHS = data(1).M2 *  V2A; % LHS
M2V2_RHS = data(1).S21*P1Aart+P2Aart; % RHS [ok]
MBV1_LHS = data(2).M1 *  V1B; % LHS
MBV1_RHS = P1Bart+data(2).S12*P2Bart;  % RHS [ok]
MBV2_LHS = data(2).M2 *  V2B; % LHS
MBV2_RHS = data(2).S21*P1Bart+P2Bart; % RHS [ok]

eps = 1e-10;
assert(meq(M1V1_LHS, M1V1_RHS, eps));
assert(meq(M2V2_LHS, M2V2_RHS, eps));
assert(meq(MBV1_LHS, MBV1_RHS, eps));
assert(meq(MBV2_LHS, MBV2_RHS, eps));
%}