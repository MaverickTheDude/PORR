function input = FI_HDCA_4links()
% Wahadlo z dwoma czlonami + sila grawitacji
L = 0.4; % [m]
mass = 0.5;
Icm = 1/12*mass*L^2;
M = diag([mass mass Icm]);
%% czlon 1
id = 1;
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
type = 'rot';
body1 = newBody(id, M, sC1, sC2, type);
%% czlon 2
id = 2;
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
type = 'rot';
body2 = newBody(id, M, sC1, sC2, type);
%% czlon 3
id = 3;
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
type = 'rot';
body3 = newBody(id, M, sC1, sC2, type);
%% czlon 4
id = 4;
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
type = 'rot';
body4 = newBody(id, M, sC1, sC2, type);

element = elementsInfo(mass, [0;0], L);
%% Ustatwienia dodatkowe
dt = 0.01;      q0 = [1; 1; 0; 0];
Tk = 2;         dq0 = [1; 2; -1; 0];

input = cell(1,7);
input{1} = [body1 body2 body3 body4];
input{3} = [element element element element];
input{7} = struct('dt',dt,'Tk',Tk,'q0',q0,'dq0',dq0,...
          'HDCAfun',@HDCAfun4R,'ACCfun',[]); 
end