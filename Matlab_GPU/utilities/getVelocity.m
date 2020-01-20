function [V1, V2] = getVelocity(body, q, dq, input)
% Returns velocity vector of n-th body (n==body) at handle H1
% (or H2) by sequential application of V2 == S'*V1 + Hq
data = getData(q, dq, input);
V1 = data(1).H*dq(1);
for j = 2:body
    V1 = data(j-1).S12'*V1 + data(j).H*dq(j);
end
V2 = data(body).S12'*V1;
end