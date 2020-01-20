function d1Kstar = recurseAcceleration(K, p, dp, q, dq, M1N, d1N, V1tab, sigmaTab, input)
% K - current body under recursion
% i == K | j == N

I = eye(3);
N = K + 1; % child of body K
M1N = M1N(:,:,N);
d1N = d1N(:,N);

data = getData(q, dq, input);
Hn = data(N).H;
dS12 = data(K).dS12;
S12k = data(K).S12;

P1Nart =  data(N).H*p(N) + data(N).D*sigmaTab(:,N);
bias = @(i) (data(i).dH * p(i) + data(i).H * dp(i) + data(N).dD*sigmaTab(:,i));
d1K = bias(K) - data(K).dM1 * V1tab(:,K) - S12k*bias(N) - dS12 * P1Nart;
coef = 1 / (Hn' * M1N * Hn);
d1Nstar = (I - coef*M1N*(Hn*Hn')) * ...
         (d1N - M1N*(dS12.'*V1tab(:,K) + data(N).dH*dq(N)));
d1Kstar =  d1K + S12k*d1Nstar;
end