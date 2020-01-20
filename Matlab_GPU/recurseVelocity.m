function [M1Kstar, b1Kstar] = recurseVelocity(p, q, M1N, b1N, K, input)
% K - current body under recursion
% i == K | j == N

I = eye(3);
N = K + 1; % child of body K
M1N = M1N(:,:,N);
b1N = b1N(:,N);

data = getS(q, input);
Hn = data(N).H;
Hk = data(K).H;
S12k = data(K).S12;
b1K = Hk*p(K) - S12k * Hn*p(N);
coef = 1 / (Hn' * M1N * Hn);
M1Nstar = M1N * S12k.' - coef * M1N*(Hn*Hn')*M1N*S12k.';
b1Nstar = (I - coef*M1N*(Hn*Hn')) * b1N;
M1Kstar =  data(K).M1 + S12k * M1Nstar;
b1Kstar =  b1K + S12k*b1Nstar;

end