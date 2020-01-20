function atTime = atTime(structure, ind)
% At time returns ind-th time instant of all sequences within structure
% Could be more general version of 'reverseTime', however now we use it to
% get state structure at final time
% used by:
% 1) sequentialAdjointSolver
% 2) AdjointHamEqs

if strcmp(ind,'end')
    ind = length(structure.T);
end

names = fieldnames(structure);
atTime = struct();

for i = 1:length(names)
    field = structure.(names{i});
    atTime = setfield(atTime, names{i}, field(:, ind));
end

