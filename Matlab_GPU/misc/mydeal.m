function varargout = mydeal(arrayIn)
% Revamp Matlabowej funkcji 'deal', zeby wygodniej przepisywac dane 
% z wektora do poszczegolnych zmiennych
% arrayIn musi szereg kulumnowych wektorow do przypisania wejsciom

if nargout ~= size(arrayIn, 2)
    error(['The number of outputs must be equal to the number '...
           'of columns of the input array'])
end

varargout = cell(nargout, 1);
for i = 1:length(varargout)
    varargout{i} = arrayIn(:,i);
end

