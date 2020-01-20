function input = FI_HDCA_2links()
% Wahadlo z dwoma czlonami + sila grawitacji
L = 0.4; % [m]
mass = 0.5;
Icm = 1/12*mass*L^2;
M = diag([mass mass Icm]);
%% czlon 1
id = 1;
type = 'rot';
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
body1 = newBody(id, M, sC1, sC2, type);
%% czlon 2
id = 2;
type = 'rot';
sC1 = [-L/2; 0]; sC2 = [L/2; 0];
body2 = newBody(id, M, sC1, sC2, type);

element = elementsInfo(mass, [0;0], L);
%% Ustatwienia dodatkowe
dt = 0.01;      q0 = [1; 1];
Tk = 2;         dq0 = [1; 2];

input = cell(1,7);
input{1} = [body1 body2];
input{3} = [element element]; 
input{7} = struct('dt',dt,'Tk',Tk,'q0',q0,'dq0',dq0,'HDCAfun',@HDCAfun2R,...
                  'ACCfun',@ACCfunCP); 
end