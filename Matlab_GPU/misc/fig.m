function fig(n)
if nargin == 1;
    if n == 1
        close all
    end
    figure(n)
else
    figure
end
set(gcf,'WindowStyle','docked','color','white');
% set(gcf,'color','white');