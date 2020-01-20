function ind = findT(t, input)
%findT returns index which corresponds to t-th time instant for a given
%input

dt = input{7}.dt;
Tk = input{7}.Tk;
T = 0:dt:Tk;
ind = find( abs(T - t) < dt/2 );

if isempty(ind)
    ind = find( abs(T - 0.75*t) < dt );
    if length(ind) > 1
        ind = pickIndex(T, t, ind);
    end
end

if  length(ind) > 1
    ind = pickIndex(T, t, ind);
end

end

function ind = pickIndex(T, t, indices)
% Wybieramy indeks najblizszy chwili czasowej t.
% Jezeli wystapi sytacja, ze length(indeces)>=3 -> trzeba uogolnic funkcje
if length(indices) >= 3
    warning('findT: indeks ind > 2')
end
eps1 = abs(T( indices(1) ) - t);
eps2 = abs(T( indices(2) ) - t);
if eps1 < eps2
    ind = indices(1);
elseif eps2 < eps1
    ind = indices(2);
elseif eps2 == eps1
    ind = indices(1);   % przyjete arbitralnie
end

end

