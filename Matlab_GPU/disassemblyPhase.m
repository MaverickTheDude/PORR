function [structA, structB] = disassemblyPhase(structC)

% currently disassembly is performed only with the revolute joints (for future reference)
D = [eye(2); [0 0]];
T1A = structC.T1;        T2B = structC.T2;
Q1Aart = structC.Q1;     Q2Bart = structC.Q2;
ksiC = structC.ksi; 
ksiA = ksiC.ksiA;        ksiB = ksiC.ksiB;
Q1Bacc = ksiC.Q1Bacc;    S12B = ksiC.S12B;

C = - D.' * (ksiB.i11 + ksiA.i22) * D; % no inverse here (!)
W = D * (C \ D.');
b = D.' * (ksiB.i10 - ksiA.i20);
beta = D * (C \ b);
T1B =  W * ksiB.i12 * T2B - W * ksiA.i21 * T1A + beta; % srodek pary kinematycznej wzory (29)
T2A =  - T1B;                                   % ze wzoru (29)

Q1Bart = Q1Bacc - S12B * Q2Bart; % wzor (76)
Q2Aart = -Q1Bart;

structA = struct('ksi',ksiA,'T1',T1A,'Q1',Q1Aart,'T2',T2A,'Q2',Q2Aart);
structB = struct('ksi',ksiB,'T1',T1B,'Q1',Q1Bart,'T2',T2B,'Q2',Q2Bart);
end

