function forcesAB = disassembleBlock(blockA, blockB, blockC, forcesC)

% currently disassembly is performed only with the revolute joints (for future reference)
D = [eye(2); [0 0]];     ksi = @(i) (3*i-2:3*i);
T1A = forcesC(1,:).';        T2B = forcesC(2,:).';
Q1Aart = forcesC(3,:).';     Q2Bart = forcesC(4,:).';

C = - D.' * (blockB(ksi(1),:) + blockA(ksi(4),:)) * D; % no inverse here (!)
W = D * (C \ D.');
b = D.' * (blockB(13,:).' - blockA(14,:).');
beta = D * (C \ b);
T1B =  W * blockB(ksi(2),:) * T2B - W  * blockA(ksi(3),:) * T1A + beta; % srodek pary kinematycznej wzory (29)
T2A =  - T1B;                                   % ze wzoru (29)

Q1Bart = blockB(end,:).' - blockB(end-3:end-1,:) * Q2Bart; % wzor (76)
Q2Aart = -Q1Bart;

forcesA = [T1A.'; T2A.'; Q1Aart.'; Q2Aart.'];
forcesB = [T1B.'; T2B.'; Q1Bart.'; Q2Bart.'];
forcesAB = [forcesA, forcesB];
end

