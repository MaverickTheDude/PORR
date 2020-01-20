function J = testFun(u, input)
% Customowa funkcja definijaca J dla wyliczenia grad(J) metoda roznic
% skonczonych. Cala zawartosc jest do zastapienia

state = DynamicHamiltonTrapz(input, u);
% addpath('old_projects\DHAM_Cart_Pendulum')
% J = cartPendCost(state, input);

T = state.T;
fiZad = pi/2;
displ = state.displ;
[n, N] = size(displ);
fi2F = displ(6,end);
S = 0.5*(fi2F - fiZad)^2;
J = S;

end

