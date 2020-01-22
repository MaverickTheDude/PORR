function Ytab = GPU_RK45(input)
% wahadla polaczone sekwencyjnie.
N = input.Nbodies;
q0 = input.q0;      p0 = input.p0;
dt = input.dt;      Tk = input.Tk;
T = 0:dt:Tk;
Nsamples = length(T);

Ytab = zeros(2*N, Nsamples);
y_m1 = [p0; q0];
Ytab(:,1) = y_m1;

for i = 2 : Nsamples
    k1 = RHS(T(i-1), y_m1, input);
	k2 = RHS(T(i-1) + dt/2, y_m1 + dt/2*k1, input);
	k3 = RHS(T(i-1) + dt/2, y_m1 + dt/2*k2, input);
	k4 = RHS(T(i-1) + dt,   y_m1 + dt*k3, input);
	y = y_m1 +  dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    Ytab(:,i) = y;
    y_m1 = y;
end
