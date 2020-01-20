function [p0, sigma0] = dq0top0(q, dq, input)
data = getData(q, dq, input);
p0 = zeros(length(q), 1);
sigma0 = zeros(2*length(q), 1);
P1B = @(p,sigma,i)(data(i).H*p + data(i).D*sigma); % P1Bart := -P2Aart
n = length(input{1});

V1 = getVelocity(n, q, dq, input);
sigma = data(n).D' * data(n).M1*V1;
sigma0(end-1:end) = sigma;
p0(n) = data(n).H' * data(n).M1*V1;
for i = n-1:-1:1
    V1 = getVelocity(i, q, dq, input);
    p0(i) = data(i).H' * (data(i).M1*V1 + data(i).S12 * P1B(p0(i+1), sigma, i+1));
    sigma = data(i).D' * (data(i).M1*V1 + data(i).S12 * P1B(p0(i+1), sigma, i+1));
    sigma0(2*i-1:2*i) = sigma;
end
end