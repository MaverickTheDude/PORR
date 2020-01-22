function [V1, V2] = getV(blockA, forcesA)
ksi = @(i) (3*i-2:3*i);

V1 = blockA(ksi(1),:)*forcesA(1,:).' + blockA(ksi(2),:)*forcesA(2,:).' + blockA(13,:).';
if nargout == 2
    V2 = blockA(ksi(3),:)*forcesA(1,:).' + blockA(ksi(4),:)*forcesA(2,:).' + blockA(14,:).';
end