function dY = RHS(t, Y, input)
% N-link pendulum in joint coordinates
% cala struktura z danymi zawarta jest w jednej macierzy:
% Assembly phase: [i11; i12; i21; i22; i10; i20; S12; Q1acc]
% Disassembly phase: [T1A; T2A, Q1Art; Q2Art]

N = input.Nbodies;
dt = input.dt;      Tk = input.Tk;
T = 0:dt:Tk;
tiers = input.tiers;
tiersInfo = input.tiersInfo;
ind = @(i) (3*i-2:3*i);

p = Y(1:N);         H = [0; 0; 1]; 
q = Y(N+1:end);     D = [eye(2); [0 0]];

AsmTree = cell(tiers);
FcsTree = cell(tiers);
for i = 1 : tiers
    AsmTree{i} = zeros(18, 3 * tiersInfo(i));
    FcsTree{i} = zeros( 4, 3 * tiersInfo(i));
end

% zamiana wsp zlaczowych na globalne
fi = cumsum(q);

for i = 1 : N-1
    AsmTree{1}(:,3*i-2:3*i) = initCoefs(fi(i), p(i), p(i+1), input);
end
AsmTree{1}(:,3*N-2:3*N) = initCoefs(fi(N), p(N), 0, input);

for i = 2 : tiers
    for j = 1 : tiersInfo(i)
        % indeksowanie po j: 2*k-1:2*k, gdzie w miejsce k wstawiamy ind(j)
        AsmTree{i} (:, ind(j)) = ...
      assembleBlocks(AsmTree{i-1}(:,6*j-5:6*j-3), AsmTree{i-1}(:,6*j-2:6*j));
    end
end


blockS = AsmTree{tiers};

% base body connection
c = - D.' * blockS(ind(1),:) * D; % no inverse here (!)
T1S = D * (c \ D.') * blockS(13,:).';     % wzor (43)
T2S = zeros(3,1);
Q1Sart = blockS(end,:).'; % articulated forces
Q2Sart = zeros(3,1);
FcsTree{tiers} = [T1S.'; T2S.'; Q1Sart.'; Q2Sart.'];

for i = tiers : -1 : 2
    for j = 1 : tiersInfo(i)
        FcsTree{i-1}(:, 6*j-5:6*j) = disassembleBlock(...
AsmTree{i-1}(:,6*j-5:6*j-3), AsmTree{i-1}(:,6*j-2:6*j),...
AsmTree{i}(:,ind(j)), FcsTree{i}(:, ind(j)));
    end
end

P1art = zeros(3,N);
P1art(:,1) = FcsTree{1}(1,ind(1)).' + H * p(1);
[dp, dq] = deal(zeros(N,1));
dq(1) = H.' * getV(AsmTree{1}(:,ind(1)), FcsTree{1}(:,ind(1)));

for i = 2 : N
    [~, V1B] = getV(AsmTree{1}(:,ind(i)), FcsTree{1}(:,ind(i)));
    V2A = getV(AsmTree{1}(:,ind(i-1)), FcsTree{1}(:,ind(i-1)));
    dq(i) = H.' * (V1B - V2A);
    P1art(:,i) = FcsTree{1}(1,ind(i)).' + H * p(i);
end
dfi = cumsum(dq);

des = zeros(3, N);
for i = N-1 : -1 : 1
    [~, dSc2] = getdS2(fi(i), dfi(i), input);
    dSc1 = getdS2(fi(i+1), dfi(i+1), input);
    des(:,i) = (dSc2 - dSc1) * P1art(:,i+1) + des(:,i+1); 
end

for i = 1 : N
    dSc1 = getdS2(fi(i), dfi(i), input);
    Q1art = FcsTree{1}(3,ind(i)).';
    dp(i) = H' * (Q1art-dSc1*P1art(:,i) + des(:,i));
end


dY = [dp; dq];


