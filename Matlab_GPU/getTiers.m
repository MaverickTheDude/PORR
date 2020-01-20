function [tiers, tiersInfo] = getTiers(Nbodies)

tiers = ceil(log2(Nbodies))+1;
tiersInfo = zeros(1,tiers);
Ntmp = Nbodies;
for i = 1 : tiers
    tiersInfo(i) = Ntmp;
    Ntmp = ceil(Ntmp / 2);
end


