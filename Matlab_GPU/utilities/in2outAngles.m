function [FI, OM] = in2outAngles(q, dq, input)
% Przelicza zlaczowe do absolutnych dla algorytmu HDCA (tylko wsprzolrzedne katowe)
calculateVelocity = true;
if nargin == 2
    calculateVelocity = false;
    input = dq;    end
[FI, OM] = deal(zeros(length(q), 1));
[fi, om] = deal(0);
for i = 1:length(input{1})
    if strcmp(input{1}(i).type,'rot')
        fi = fi + q(i);
        if calculateVelocity, om = om + dq(i); end
    end
    FI(i) = fi;
    OM(i) = om;
    if nargout == 1 % safety measure
        OM = [];    end
end