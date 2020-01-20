function [ped, dped] = getPhysicalP(state, input)
% getPhysicalP oblicza fizyczne ped = M dq i jego pochodna dp
% dla sformulowania klasycznego (globalnego)

[n, N] = size(state.displ);
[pedT, ped, dped] = deal(zeros(n, N));

for i = 1 : N
    q = state.displ(:,i);
    dq = state.vel(:,i);
    d2q = state.acc(:,i);
    Fq = Jacobian(q, input);
    M1 = Mass(q, input);
    dM1 = MassDot(q, dq, input);
    pedT(:,i) = state.ped(:,i) - Fq.'*state.sigma(:,i);
    ped(:,i) = M1 * dq;
    dped(:,i) = dM1*dq + M1*d2q;
    assert( meq(pedT(:,i),ped(:,i)) )
end