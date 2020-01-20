addpath('inputs','utilities','misc')
input = FI_HDCA_4links();
q0 = [0; -pi/4; 0; pi/4]; dq0 = [-1.0; 2; 0; 0.5]; % 4 dof
dt = 0.01;          Tk = 2;
input{7}.q0 = q0;     input{7}.dq0 = dq0;
input{7}.dt = dt;     input{7}.Tk = Tk;
stateHDCA = HDCASolver(input);
stateH = joint2absolute(stateHDCA, input);
calculateEnergy(stateH, input);
%% wyliczenie p0
p0 = dq0top0(input{7}.q0, input{7}.dq0, input);

%% nasz main
Nbodies = 4;
[tiers, tiersInfo] = getTiers(Nbodies);
L = 0.4; % [m]
mass = 0.5;
Icm = 1/12*mass*L^2;
M = diag([mass mass Icm]);
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
input = struct('Nbodies',Nbodies,'q0',q0,'p0',p0,'dt',dt,'Tk',Tk,...
    'tiers',tiers,'tiersInfo',tiersInfo,'M',M,'sC1',sC1,'sC2',sC2);

%%
q0 = [0; -pi/4; 0; pi/4];

fi = arrayfun(@cumsum,q0)
