function dS = getdS(q, dq, input)
% returns auxaliary matrices and vector for current state based on the input
% shift Matrix - dwie redefinicje:
% 1) S := S' (zeby byc "on the same page" z artykulem Pawla)
% 2) squeezed_tilde := -squeezed_tilde:
%    LHS: wektor s jest od biezacego c-sys do innego
%    RHS: wektor s jest od jakiegos c-sys (srodka masy) do biezacego
if length(q) ~= length(input{1})
    error('funkcja getdS przyjmuje tylko wektor wsp ZLACZOWYCH');
end
I = eye(3); Om = [0, -1; 1, 0]; 
[fi, omega] = in2outAngles(q, dq, input);

dS = struct('dS12',[],'dS1c',[],'dS2c',[],...
              'dS21',[],'dSc1',[],'dSc2',[],...
              'dH',[],'dD',[],'om',[]);
N = length(input{1});
for i = 1 : N
    [dS12, dS21, dS1c, dSc1, dS2c, dSc2] = deal(zeros(3));
    sC1_loc = input{1}(i).sC1;
    if i < N && strcmp(input{1}(i+1).type,'trans')
        % zakladamy, ze wsp. niezalezna wyznacza kierunek zgodny z osia Ox
        % poprzedniego (czyli obecnego, i-tego) c-sys
        s12_loc = [q(i+1); 0];
        sC2_loc = sC1_loc + s12_loc;
    else
        s12_loc = input{1}(i).s12;
        sC2_loc = input{1}(i).sC2;
    end

    s12 = Rot(fi(i))*s12_loc;   s21 = -Rot(fi(i))*s12_loc;
    sC1 = Rot(fi(i))*sC1_loc;   s1C = -Rot(fi(i))*sC1_loc;
    sC2 = Rot(fi(i))*sC2_loc;   s2C = -Rot(fi(i))*sC2_loc;
    
    om = omega(i);
    dS12(3,1:2) = (-om*s12).';   dS21(3,1:2) = (-om*s21).';
    dSc1(3,1:2) = (-om*sC1).';   dS1c(3,1:2) = (-om*s1C).';   
    dSc2(3,1:2) = (-om*sC2).';   dS2c(3,1:2) = (-om*s2C).';   
    % poprawka tylko w strone 1 -> 2, c -> 2 wynika z zalozenia,
    % ze mamy pare wahadlo->suwak (nic innego nie przewidujemy).
    % Sc1 i dSc1 sa bez zmian, bo suwak jest z prawej strony.
    if i < N && strcmp(input{1}(i+1).type,'trans')
        ds = [dq(i+1); 0];
        dS12(3,1:2) = dS12(3,1:2) + (Om*Rot(fi(i))*ds).';
        dSc2(3,1:2) = dSc1(3,1:2) + dS12(3,1:2);
    end
    
    if strcmp(input{1}(i).type, 'trans')
        u = Rot(fi(i))*[1;0];
        dH = [Om*u*om; 0];
        dD = [[-u*om;0], zeros(3,1)];
    else
        dH = zeros(3,1);    dD = zeros(3,2);
    end
    
    dS(i) = struct('dS12',dS12,'dS1c',dS1c,'dS2c',dS2c,...
                   'dS21',dS21,'dSc1',dSc1,'dSc2',dSc2,...
                   'dH',dH,'dD',dD,'om',omega(i));
end
end