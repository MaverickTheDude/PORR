function [Em, Ek, Ep, W] = calculateEnergy(state, input)
% Oblcza calkowita energie mechaniczna oraz prace wykonana przez zewnetrzne
% sily dla sformulowania we wspolrzednych absolutnych.
% Updated: uwzgledniono uklady w ktorych c-sys nie pokrywa sie z CM
Q = state.displ;    g = 9.80665;
DQ = state.vel;     Om = Rot(pi/2);
[n, N] = size(Q);
[Ek, Ep, W] = deal(zeros(1, N));
for i = 1:N
    for j = 1:n/3
        s1C = -input{3}(j).sCM;
        q = Q(3*j-2:3*j, i);    m = input{3}(j).mass;
        dq = DQ(3*j-2:3*j, i);  J = input{3}(j).inertia;
        rCM = q(1:2) + Rot(q(3))*s1C;
        drCM = dq(1:2) + Om*Rot(q(3))*s1C*dq(3);
        Ek(i) = Ek(i) + 0.5 * m * drCM.'*drCM + 0.5 * J * dq(3)^2;
        Ep(i) = Ep(i) + m * g * rCM(2); % to do: uwzglednic potencjal sprezyny
    end
    % Work of applied forces
    if isempty(input{6})
        continue;
    end;
    for j = 1:length(input{6})
        if input{6}(j).control
            F = input{6}(j).F(i);
        else
            t = (i-1)*input{7}.dt;  F = input{6}(j).F(t);
        end
        s1C = -input{3}(j).sCM;
        id = input{6}(j).elements(1);
        fi = state.displ(3*id, i);
        r = state.displ(3*id-2:3*id-1, i) + Rot(fi)*s1C;
        u = input{6}(j).u;
        if strcmp(u, 'moment')
            M = F;
            W(i) = W(i) + M * fi;
        else
            sC1 = input{3}(id).sCM;
            F0 = F * Rot(fi) * u; % Force vector in global c-sys
            M = (Om*sC1).'*F0;
            W(i) = W(i) + F0.'*r + M*fi;
        end
    end
end

Em = Ek + Ep - W;

figure('WindowStyle','docked','color','white')
if any(W)
    plot(state.T, Ek, state.T, Ep, state.T, Em, state.T, W);
    legend('E_k', 'E_p', 'E_m', 'W','Location','best','orientation','horizontal')
else
    plot(state.T, Ek, state.T, Ep, state.T, Em);
    legend('E_k', 'E_p', 'E_m','Location','best','orientation','horizontal')
end
title('energy in the system')
ylabel energy; grid on