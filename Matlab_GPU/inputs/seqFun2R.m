function [dY, dq, dp, sigma] = seqFun2R(t, Y, input)
% Wahadlo o dwoch stopniach swobody w podejsciu sekwencyjnym.
p = Y(1:2);     H = [0; 0; 1];      I = eye(3);
q = Y(3:4);     D = [eye(2); [0 0]];
data = getS(q, input);
b1B = H * p(2);
b1A = H*p(1) - data(1).S12 * H*p(2);
coef = 1 / (H' * data(2).M1 * H);
M1Bstar = data(2).M1*data(1).S12.' - coef * data(2).M1*(H*H')*data(2).M1*data(1).S12.';
b1Bstar = (I - coef*data(2).M1*(H*H')) * b1B;
M1Astar =  data(1).M1 + data(1).S12 * M1Bstar;
b1Astar =  b1A + data(1).S12*b1Bstar;

dq1 = (H'*M1Astar*H) \ (H'*b1Astar);
V1A = H * dq1;
dq2 = (H'*data(2).M1*H) \ (H'*b1B - H'*data(2).M1*data(1).S12.'*V1A);

V1B = data(1).S12'*V1A + H*dq2;
T1A = M1Astar*V1A - b1Astar;
T1B = data(2).M1*V1B - b1B;
sigma = [T1A(1:2); T1B(1:2)];
% LHS = data(1).M1*V1A;
% RHS = T1A - data(1).S12*T1B + b1A;

dq = [dq1; dq2];

P1Aart = T1A + H*p(1);
P1Bart = T1B + H*p(2);

dS = getdS(q, dq, input);
Q1 = setForcesAtH1(q, input);
[Q1A, Q1B] = mydeal(Q1);
Q1C = Q1A + data(1).S12*Q1B;
dSC1_A = dS(1).dSc1;   dSC2_A = dS(1).dSc2;
dSC1_B = dS(2).dSc1;   dSC2_B = dS(2).dSc2;
des1 = H'*(dSC2_A - dSC1_B)*P1Bart;

dp1 = H' * (Q1C - dSC1_A*P1Aart) + des1;
dp2 = H' * (Q1B - dSC1_B*P1Bart);
dp = [dp1; dp2];

dY = [dp; dq];

% Sformulowanie accumulated (zamiast articulated) dla dp1
% P1A = data(1).M1 * V1A;
% P1B = data(2).M1 * V1B; 
% Q1A = Q1(:,1);
% Hp2 = Q1Bart - dSC1_B * P1Bart;
% dp1 = H' * (Q1A - dSC1_A*P1A + dS(1).dS12*P1B + data(1).S12* Hp2 );