function input = FI_HDCA_CP(u, q0, dq0)
% Wahadlo na wozku + sila grawitacji (sformulowanie HDCA / Seq)
if nargin == 1
    q0 = [0; -pi/2];
    dq0 = [0; 0];
end

L = 1; h = L/2; % [m]
mC = 0.8; mP = 0.4;
Icmc =  1/12*mC*(L^2 + h^2); Icmp = 1/12*mP*L^2;
Mc = diag([mC mC Icmc]); Mp = diag([mP mP Icmp]);
%% czlon 1
id = 1;
sC1 = [0; 0];      sC2 = [0; 0];
type = 'trans';    v = [0; 1];
s1P = [1; 0];      
% F = [forceVal; 0]; T = 0;
Z = zeros(1,length(u));
F = [u.'; Z]; T = Z;
forceControl = struct('type','control','s1P',s1P,'F',F,'T',T);
c = 0.5;
% forceDamp = struct('type','dampingTrans','s1P',[],'F',c,'T',[]);
forceDamp = [];
force = [forceControl, forceDamp];
body1 = newBody(id, Mc, sC1, sC2, type, v, force);
%% czlon 2
id = 2;
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
type = 'rot';
c = 0.2;
% force = struct('type','dampingRot','s1P',[],'F',c,'T',[]);
force = [];
body2 = newBody(id, Mp, sC1, sC2, type, [], force);
%{
%% czlon 3
id = 3; type = 'rot';
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
body3 = newBody(id, Mp, sC1, sC2, type);
%% czlon 4
id = 4; type = 'rot';
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
body4 = newBody(id, Mp, sC1, sC2, type);
%}
%% Si³a zadana
element = 1; u_vec = [1; 0];
F = u;
sAi = [0; 0]; control = true;
force = fillForce(element, u_vec, F, sAi, control); 
%{
element = 1; u_vec = [1; 0];
F = @(t)(forceVal);
sAi = [0; 0]; control = false;
force = fillForce(element, u_vec, F, sAi, control); 
%}
%% Czlony
Cart = elementsInfo(mC, [0;0], L, h);
Pend = elementsInfo(mP, [0;0], L);
%% Ustatwienia dodatkowe
dt = 0.02;
Tk = 0.5;
% q0 = [0; -pi/2; 0; 0];  dq0 = [0; 0; 0; 0];
% q0 = [0; -pi/2];          dq0 = [0; 0];

input = cell(1,7);
% input{1} = [body1 body2 body3 body4];
% input{3} = [Cart Pend Pend Pend];
input{1} = [body1 body2];
input{3} = [Cart Pend];
input{6} = force;
input{7} = struct('dt',dt,'Tk',Tk,'q0',q0,'dq0',dq0,'HDCAfun',@HDCAfunCP,...
                  'ACCfun',@ACCfunCP,'dS',@dS);  % @ACCfunCP
end


function [blank, ETA_KSI_MI0] = dS(endState, input)
blank = [];
q = endState.displ; p = endState.ped;
v = endState.vel;   Om = [0, -1; 1, 0];
I = eye(2);         
n = length(q);      
fiZad = pi/2;       fi2 = q(6);
Sq = [zeros(1,3), 0, 0, fi2-fiZad];
% Sv = zeros(1,n);
Sp = zeros(1,n);
vZad = 0;
Sv = [zeros(1,3), 0, 0, v(6)-vZad];
M = Mass(q, input);
Fq = Jacobian(q, input);
m = 2*length(input{1}) + 2*length(input{2});
n = length(q); z = zeros(1,n);
ksi0 = -Sp.';

A = Q_v(q, input) + DFqTSIG_DQ(q, endState.sigma, input); % bez czlonu zwiazanego z przesunieciem CM
DFq = DF_q(q, endState.vel, input);
eta0_mi0 = [M, Fq.'; Fq, zeros(m)] \ [Sq.' - A.'*ksi0; -DFq * ksi0];
ETA_KSI_MI0 = [eta0_mi0(1:n); ksi0; eta0_mi0(n+1:end)];
end