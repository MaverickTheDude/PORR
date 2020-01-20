function [state, ind] = reverseTime(tau, stateIn, input)
% Funkcja do znajdowania chwili tau w dziedzinie T przy rozwiazywania
% rownan adjointowych 'wstecz'. Wersja zewnetrzna, ale na razie
% wykorzystywana przez:
% 1) AbsAdj.HDCA_AdjointSolver (absolute coordinates).
% 2) sequentialAdjointSolver
% #1 Usunieto dLambda i d2Lambda.
T = stateIn.T;
dt = input{7}.dt; % krok czasowy wzglêdem zadania dynamiki, ¿eby nie by³o kilku indeksów
Tk = input{7}.Tk;
t = Tk - tau;

ind = find( abs(T - t) < dt/2 );
if length(ind) > 1
    warning('indeks ind > 1')
    ind2 = find( abs(T - t) < dt/4 );
    if length(ind2) == 1
        ind = ind2;
    else
        ind = ind(1);
    end
end
if isempty(ind)
    ind = find( abs(T - t) < dt );
    if length(ind) > 1
        ind = ind(1);
    end
end
if isempty(ind)
    ind = find( abs(T - t) < 2*dt );
    if length(ind) > 1
        ind = ind(1);
    end
end

try
    state.T = T(ind);
    state.q = stateIn.displ(:,ind);      state.dq = stateIn.vel(:,ind);
    state.sigma = stateIn.sigma(:,ind);  state.lambda = stateIn.lambda(:,ind);
    state.ped = stateIn.ped(:,ind);      %state.dped = stateIn.dped(:,ind); -> niepotrzebne w SeqAdj
    state.d2q = stateIn.acc(:,ind);     
catch ME 
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
    state.T = T(ind);                
    state.q = stateIn.Q(:,ind);          state.dq = stateIn.dQ(:,ind);
    state.sigma = stateIn.sigma(:,ind);  state.lambda = stateIn.lambda(:,ind);
    state.ped = stateIn.pJoint(:,ind);     state.dped = stateIn.dpJoint(:,ind);
    state.eps = stateIn.eps(:,ind);     
    else
        error('cos nie tak')
    end
end

