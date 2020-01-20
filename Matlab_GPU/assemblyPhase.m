function [Q1Cacc, ksiC, S12C] = assemblyPhase(assemblyA, assemblyB)
% Assembles the force C1 acting on first handle of compound body C
% Calculates coefficients of body C for given coefs of bodies A & B

S12C = assemblyA.S12*assemblyB.S12;
% accumulated force in handle 1 of body C
Q1Cacc = assemblyA.Q1 + assemblyA.S12 * assemblyB.Q1; % LHS(73) = LHS(74)
ksiA = assemblyA.ksi;   ksiB = assemblyB.ksi;

% currently assembly is performed only with the revolute joints (for future reference)
D = [eye(2); [0 0]]; % how to generelize this?

C = - D.' * (ksiB.i11 + ksiA.i22) * D; % no inverse here (!)
W = D * (C \ D.');
b = D.' * (ksiB.i10 - ksiA.i20);
beta = D * (C \ b);
ksi11 = ksiA.i11 + ksiA.i12 * W * ksiA.i21;
ksi22 = ksiB.i22 + ksiB.i21 * W * ksiB.i12;
ksi12 = -ksiA.i12 * W * ksiB.i12;
ksi21 = -ksiB.i21 * W * ksiA.i21;
ksi10 = ksiA.i10 - ksiA.i12 * beta;
ksi20 = ksiB.i20 + ksiB.i21 * beta;
ksiC = struct('i11',ksi11,'i22',ksi22,'i12',ksi12,...
              'i21',ksi21,'i10',ksi10,'i20',ksi20,...
 'ksiA',ksiA,'ksiB',ksiB,'S12B',assemblyB.S12,'Q1Bacc',assemblyB.Q1);
assert(all(all(abs(ksi21 - ksi12.') < 1e-10))) % security assertion
if nargout == 1     % overload dla jednego wyjscia
    Q1Cacc = struct('Q1',Q1Cacc,'ksi',ksiC,'S12',S12C);
end