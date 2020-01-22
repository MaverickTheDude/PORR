% custom Runge Kutta w c++ vs w Matlabie
PATH = 'D:\Code\PORR\Sekwencyjny';
working_dir = cd(PATH);
file = 'results.txt';
data = importdata(file);
Nbodies = (size(data,1)-1)/2;

T = data(1,:);
pTab = data(2:Nbodies+1,:);
qTab = data(end-Nbodies+1:end,:);

cd(working_dir);

%% Matlab
% q0 = [pi/2; -pi/4]; dq0 = [-2; 4];
% input = FI_HDCA_2links(); 
q0 = [0; -pi/4; 0; pi/4]; dq0 = [-1.0; 2; 0; 0.5]; % 4 dof
dt = 0.01;          Tk = 1;
input = FI_HDCA_4links();
input{7}.q0 = q0;     input{7}.dq0 = dq0;
input{7}.dt = dt;     input{7}.Tk = Tk;
stateHDCA = HDCASolver(input);
stateH = joint2absolute(stateHDCA, input);
%% wyniki
b = 2;
fig; plot(T, qTab(b,:)); hold on;
plot(stateHDCA.T, stateHDCA.Q(b,:));
legend('c++','Matlab'); grid on; 

fig; plot(T, pTab(b,:)); hold on;
plot(stateHDCA.T, stateHDCA.pJoint(b,:));
legend('c++','Matlab'); grid on;
%% Wykres
fig;
for i = 1:4
    subplot(2,2,i)
    plot(T, qTab(i,:), stateHDCA.T, stateHDCA.Q(i,:));
    ylabel(['$\varphi_', num2str(i), '$'],'interpreter','latex');
    xlabel t
    legend('c++','Matlab','Location','Best'); grid on; 
end