addpath('inputs','utilities','misc')
%% nasz main
Nbodies = 2048;
[tiers, tiersInfo] = getTiers(Nbodies);
L = 0.04; % [m]
mass = 0.05;
Icm = 1/12*mass*L^2;
M = diag([mass mass Icm]);
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
q0 = zeros(Nbodies,1); q0(1) = pi/6;
p0 = zeros(Nbodies,1); 
dt = 0.01; Tk = dt;
% p0 = dq0top0(input{7}.q0, input{7}.dq0, input);
inputTMP = struct('Nbodies',Nbodies,'q0',q0,'p0',p0,'dt',dt,'Tk',Tk,...
    'tiers',tiers,'tiersInfo',tiersInfo,'M',M,'sC1_loc',sC1,'sC2_loc',sC2);

T = 0:inputTMP.dt:inputTMP.Tk;
tic
Ytab = GPU_RK45(inputTMP);
t = toc
pTab = Ytab(1:Nbodies,:);
qTab = Ytab(Nbodies+1:end,:);

%% compare
input = FI_HDCA_4links();
q0 = [0; -pi/4; 0; pi/4]; dq0 = [-1.0; 2; 0; 0.5]; % 4 dof
dt = 0.01;          Tk = 2;
input{7}.q0 = q0;     input{7}.dq0 = dq0;
input{7}.dt = dt;     input{7}.Tk = Tk;
stateHDCA = HDCASolver(input);
stateH = joint2absolute(stateHDCA, input);
calculateEnergy(stateH, input);

%% wykresy
x = [512 1024 1500 2048]
y = [0.911 1.807 2.7 3.586]

fig; plot(x,y,'o'); grid on
