function constraintPosition = getContraintPosition(i, j, inputC)
% getContraintPosition returns the constraint's id (position in Phi vector)
% for a pair of bodied with ids 'i' and 'j'

% rotational joints
for ii = 1 : length(inputC{1})
    elements = sort(inputC{1}(ii).elements);
    if  all(elements == sort([i, j]))
        constraintPosition = inputC{1}(ii).ConstrPosition;
        return
    end
end
% translational joints
for ii = 1 : length(inputC{2})
    elements = sort(inputC{2}(ii).elements);
    if  all(elements == sort([i, j]))
        constraintPosition = inputC{2}(ii).ConstrPosition;
        return
    end
end
% no common joints found - throw an error
error('no common joints found')
end