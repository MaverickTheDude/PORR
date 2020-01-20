function state = joint2absolute(stateH, input)
% Transfers joint coordinates into abolsute ones for HDCA projects.
% currently: Translational joints need to be taken into account separately

dQ = stateH.dQ;     Tk = input{7}.Tk;
Q = stateH.Q;       dt = input{7}.dt;
T = 0:dt:Tk;        state.sigma = stateH.sigma;
state.T = T;
if ~isempty(input{7}.ACCfun)
    EPS = stateH.eps; end

[b, N] = size(Q);
[displ, vel, ped, acc] = deal(zeros(3*b, N));

for i = 1:N
    q = Q(:, i);
    dq = dQ(:, i);
    if ~isempty(input{7}.ACCfun)
        eps = EPS(:,i); 
        d2q = out2inAngles(eps, input);
    end
    data = getData(q, dq, input);
    for body = 1:b
        V1 = getVelocity(body, q, dq, input);   S1C = data(body).S1c;
        M1 = data(body).M1;                     SC1 = data(body).Sc1;
        P1 = M1 * V1;
        vel(3*body-2:3*body, i) = S1C.' * V1;
        ped(3*body-2:3*body, i) = SC1   * P1;
        displ(3*body-2:3*body, i) = getPosition(body, q, input);
        if ~isempty(input{7}.ACCfun)
            acc(3*body-2:3*body, i) = ...
            getAcceleration(body, q, dq, d2q, input);
        end
    end
end

state.displ = displ;
state.vel = vel;
state.ped = ped;
state.acc = acc;
if isfield(stateH, 'lambda')
    state.lambda = stateH.lambda; end
% if ~isempty(input{7}.ACCfun)
%     state.acc = []; end
end

function d2q = out2inAngles(eps, input)
% Zamiana przyspieszen zewnetrznych do wewnetrznych. Powod: obliczajac
% przyspieszenia ze wzoru d2q = H' * dV otrzymujemy dV(3) jako wspolrzedna ABSOLUTNA. 
% Wersja customowa dla multi-link wahadla, albo poziomego wozka + wahadla
d2q = diff([0; eps]);
if strcmp(input{1}(1).type,'trans')
    d2q = [eps(1); diff([0; eps(2:end)])];  end
end

function [dVc, dV1, dV2] = getAcceleration(body, q, dq, d2q, input)
% Returns acceleration vector of n-th body (n==body) at CM, handle H1
% (or H2) by sequential application of d/dt V2 == S'*V1 + Hq
data = getData(q, dq, input);
V1 = data(1).H*dq(1);
dV1 = data(1).H*d2q(1);
for j = 2:body
    dV1 = data(j-1).S12.'*dV1 + data(j-1).dS12.'*V1 + data(j).H*d2q(j);
    V1 = data(j-1).S12'*V1 + data(j).H*dq(j);
end
dVc = data(body).S1c.'*dV1 + data(body).dS1c.'*V1;
dV2 = data(body).S12.'*dV1 + data(body).dS12.'*V1;
end

function qAbs = getPosition(body, q, input)
% Przelicza wektor wsporzednych zlaczowych na wektor wspolrzednych
% absolutnych dla czlonu "body"
% obecnie: para postepowa tylko w wersji poziomego suwaka na poczatku lancucha
data = getS(q, input);
fi = in2outAngles(q, input);
r = zeros(2,1);
if strcmp(input{1}(1).type,'trans')
    r = [q(1); 0];           end
for i = 1:body-1
    r = r + data(i).s12;
end
r = r + data(body).s1C;
qAbs = [r; fi(body)];
end