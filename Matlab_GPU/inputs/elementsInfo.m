function info = elementsInfo(m,sCM,length,width)
% Nowa wersja: przetrzymuje masê, bezw³adnoœæ i sCM (z CM do O)
% Oblicza inercjê po podstawie d³ugoœæi i szerokoœci (domyœlnie 0)
if nargin == 3
    width = 0;
end
I = m*(length^2+width^2)/12;
info = struct('mass',m,'inertia',I,'sCM',sCM);
end