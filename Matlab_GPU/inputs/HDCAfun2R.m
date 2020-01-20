function [dY, dq, dp, sigma] = HDCAfun2R(t, Y, input)
% Wahadlo o dwoch stopniach swobody
p = Y(1:2);     H = [0; 0; 1];
q = Y(3:4);     D = [eye(2); [0 0]];
data = getS(q, input);

% assembly
ksi = ksiCoefs(q, p, input);
Q1 = setForcesAtH1(q, input);
assemblyA = struct('Q1',Q1(:,1),'ksi',ksi(1),'S12',data(1).S12);
assemblyB = struct('Q1',Q1(:,2),'ksi',ksi(2),'S12',data(2).S12);
[Q1Sacc, ksiS, S12S] = assemblyPhase(assemblyA, assemblyB);

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
[structA, structB] = disassemblyPhase(structS);

Q1Aart = structA.Q1;   T1A = structA.T1;   T2A = structA.T2;
Q1Bart = structB.Q1;   T1B = structB.T1;   T2B = structB.T2;
sigma = [T1A(1:2); T1B(1:2)];

% Velocity calculation
% wzor (46) zgadza sie ze wzorami (15), (16).
[V1A, V2A] = calculateVelocities(ksi(1), T1A, T2A); % wzory (15), (16)
[V1B, V2B] = calculateVelocities(ksi(2), T1B, T2B);

dq = zeros(2,1);
dq(1) = H' * V1A;
dq(2) = H' * (V1B - V2A);

% Articulated momenta
P1Aart = T1A + H*p(1);
P1Bart = T1B + H*p(2);
P2Aart = T2A - H*p(2);

dS = getdS(q, dq, input);
dSC1_A = dS(1).dSc1;   dSC2_A = dS(1).dSc2;
dSC1_B = dS(2).dSc1;   dSC2_B = dS(2).dSc2;
des1 = H'*(dSC2_A - dSC1_B)*P1Bart;

dp = zeros(2,1);
dp(1) = H' * (Q1Aart - dSC1_A * P1Aart) + des1;
dp(2) = H' * (Q1Bart - dSC1_B * P1Bart);
dY = [dp; dq];

% Sformulowanie accumulated (zamiast articulated) dla dp1
P1A = data(1).M1 * V1A;
P1B = data(2).M1 * V1B; 
Q1A = Q1(:,1);
Hp2 = Q1Bart - dSC1_B * P1Bart;
dp1 = H' * (Q1A - dSC1_A*P1A + dS(1).dS12*P1B + data(1).S12* Hp2 );

