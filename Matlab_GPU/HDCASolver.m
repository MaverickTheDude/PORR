function stateHDCA = HDCASolver(input)
% HDCA Solver.
% Obecnie: non-general purpose. Customowa funkcje definiujemy w oddzielnym
% pliku a wskaznik na nia przekazujemy w input{7}.
HDCAfun = input{7}.HDCAfun;
q0 = input{7}.q0;       dt = input{7}.dt;
dq0 = input{7}.dq0;     Tk = input{7}.Tk;
p0 = dq0top0(q0, dq0, input);
options = odeset('RelTol', 1e-7, 'AbsTol' , 1e-7);
[T, solution] = ode45(@(t,Y)HDCAfun(t,Y,input), 0:dt:Tk, [p0; q0], options);
[N, n] = size(solution);    n = n/2;
ptab = solution(:,1:n)';    Q = solution(:,n+1:end)';
dQ = zeros(n, N);           sigma = zeros(2*n, N); % kazdy czlon dodaje 1 DOF, stad uproszczenie ze m:=2*n
dptab = zeros(n, N);        
for i = 1:length(T)
    [~, dq, dp, sig] = HDCAfun(T(i), solution(i,:)', input);
    dQ(:,i) = dq;
    sigma(:,i) = sig;
    dptab(:,i) = dp;
end

[~, lambda] = diff24(sigma, T.');

stateHDCA.lambda = lambda;
stateHDCA.Q = Q;
stateHDCA.dQ = dQ;
stateHDCA.pJoint = ptab;
stateHDCA.sigma = sigma;
stateHDCA.dpJoint = dptab;
stateHDCA.T = T.';

if ~isempty(input{7}.ACCfun)
    eps = zeros(n,N);       % note: "eps(i)" to w istocie dVc(3*i)
    dVc = zeros(6,N);       % zostawiamy na wszelki wypadek
    ACCfun = input{7}.ACCfun;
    for i = 1:N
        [eps(:,i), dVc(:,i)] = ACCfun(i, stateHDCA, input);
    end
    stateHDCA.eps = eps;
    stateHDCA.dVc = dVc;
end

end