function body = newBody(id, M, sC1, sC2, type, v, force)
% gather local info about body for HDCA. All info in local c-sys fixed at the CM
% H and D describe the kinematic pair which body 'id' makes with PREVIOUS member
% currently: v given in global c-sys and constant (the motion axis is contant as well)
if strcmp(type,'rot')
    H = [0; 0; 1]; D = [eye(2); [0, 0]];
    v = [];
elseif strcmp(type,'trans')
    H = [1; 0; 0]; D = [ [v;0], [0;0;1] ];
end
if nargin < 7
    force = []; end

body = struct('id',id, 'M',M, 'sC1',sC1, 'sC2',sC2, 's12',-sC1+sC2,...
               'type',type, 'v',v, 'H',H, 'D',D, 'force',force);
end