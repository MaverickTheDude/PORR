function blockC = assembleBlocks(blockA, blockB)

S12C = blockA(end-3:end-1,:)*blockB(end-3:end-1,:);
% accumulated force in handle 1 of body C
Q1C = blockA(end,:).' + blockA(end-3:end-1,:) * blockB(end,:).'; % LHS(73) = LHS(74)

% currently assembly is performed only with the revolute joints (for future reference)
D = [eye(2); [0 0]]; % how to generelize this?
ksi = @(i) (3*i-2:3*i);

C = - D.' * (blockB(ksi(1),:) + blockA(ksi(4),:)) * D; % no inverse here (!)
W = D * (C \ D.');
b = D.' * (blockB(13,:) - blockA(14,:)).';
beta = D * (C \ b);
ksi11 = blockA(ksi(1),:) + blockA(ksi(2),:) * W * blockA(ksi(3),:);
ksi22 = blockB(ksi(4),:) + blockB(ksi(3),:) * W * blockB(ksi(2),:);
ksi12 = -blockA(ksi(2),:) * W * blockB(ksi(2),:);
ksi21 = -blockB(ksi(3),:) * W * blockA(ksi(3),:);
ksi10 = blockA(13,:).' - blockA(ksi(2),:) * beta;
ksi20 = blockB(14,:).' + blockB(ksi(3),:) * beta;

assert(all(all(abs(ksi21 - ksi12.') < 1e-10))) % security assertion

blockC = [ksi11; ksi12; ksi21; ksi22; ksi10.'; ksi20.'; S12C; Q1C.'];