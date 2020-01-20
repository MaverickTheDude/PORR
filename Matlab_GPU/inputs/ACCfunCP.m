function [d2q, dVc] = ACCfunCP(j, stateH, input)
% Funkcja oblicza przyspieszenia zlaczowe w H1 oraz H2 w zadaniach:
% wozek + wahadlo oraz 2-link wahadlo
q  = stateH.Q(:,j);         sigma = stateH.sigma(:,j);
dq = stateH.dQ(:,j);        lambda = stateH.lambda(:,j);
p  = stateH.pJoint(:,j);    ind = @(i)(2*i-1:2*i);
dp = stateH.dpJoint(:,j);

data = getData(q, dq, input);
Q1tab = setForcesAtH1(q, input, j);
i = 2; % i-th body
Q1 = Q1tab(:,i);
V1 = getVelocity(i, q, dq, input);
P1 = data(i).H*p(i) + data(i).D*sigma(ind(i));
dP1 = data(i).H*dp(i) + data(i).D*lambda(ind(i));
MdV1B = - data(i).dM1*V1 + dP1;
d2qB = data(i).H.' * (data(i).M1 \ MdV1B);

dV1B = data(i).M1 \ MdV1B;
dVcB = data(i).dS1c.'*V1 + data(i).S1c.'*dV1B;
% sformulowanie, ktore nie dziala
MV1B_wrong = - data(i).dM1*V1 + Q1 - data(i).dSc1*P1 + data(i).D*lambda(ind(i));
d2qB_wrong = data(i).H.' * (data(i).M1 \ MV1B_wrong);

i = 1;
Q1 = Q1tab(:,i);
V1 = getVelocity(i, q, dq, input);
P1B = P1;
dP1B = dP1;
P1 = data(i).H*p(i) + data(i).D*sigma(ind(i));
dP1 = data(i).H*dp(i) + data(i).D*lambda(ind(i));
MdV1A = - data(i).dM1*V1 + dP1 - data(i).dS12*P1B - data(i).S12*dP1B;
d2qA = data(i).H.' * (data(i).M1 \ MdV1A);

dV1A = data(i).M1 \ MdV1A;
dVcA = data(i).dS1c.'*V1 + data(i).S1c.'*dV1A;
% sformulowanie, ktore nie dziala
lambdas = data(i).D*lambda(ind(i)) - data(i).S12*data(i+1).D*lambda(ind(i+1));
MV1A_wrong = - data(i).dM1*V1 + Q1 - data(i).dSc1*P1 + lambdas - data(i).dS12*P1B;
d2qA_wrong = data(i).H.' * (data(i).M1 \ MV1A_wrong);

dVc = [dVcA; dVcB];
d2q = [d2qA; d2qB];
end

