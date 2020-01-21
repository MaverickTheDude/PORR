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
    'tiers',tiers,'tiersInfo',tiersInfo,'M',M,'sC1_loc',sC1,'sC2_loc',sC2);
%% sequential ex
fi = pi/2;
L = 0.4;                mass = 0.5;
p = 0.5;                pnext = 0.1;
sC1 = [-L/2; 0];        sC2 = [L/2; 0];
Icm = 1/12*mass*L^2;    M = diag([mass mass Icm]);
input = struct('M',M,'sC1_loc',sC1,'sC2_loc',sC2);
block = createBlock(fi,p, pnext,input)
%%
s = 1000;
GlobalAsm = zeros(17, 3 * s);
tic
for i = 1 : s
    GlobalAsm(:,3*i-2:3*i) = createBlock2(0.5*i,i/2, i,input);
end
t2 = toc
%% stack overflow
fi = pi/2;
L = 0.4;                mass = 0.5;
p = 0.5;                pnext = 0.1;
sC1 = [-L/2; 0];        sC2 = [L/2; 0];
Icm = 1/12*mass*L^2;    M = diag([mass mass Icm]);
M = gpuArray(M);
sC1 = gpuArray(sC1);
sC2 = gpuArray(sC2);
input = struct('M',M,'sC1_loc',sC1,'sC2_loc',sC2);
block = createBlock(fi,p, pnext,input)


%%
s = 1000;
GlobalAsm = zeros(17, 3 * s, 'gpuArray');
tic
for i = 1 : s
    i = gpuArray(i);
    GlobalAsm(:,3*i-2:3*i) = createBlock(0.5*i,i/2, i,input);
end
t = toc