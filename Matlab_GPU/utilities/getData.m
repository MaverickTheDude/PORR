function data = getData(q, dq, input)
% returns auxaliary matrices and vector for current state based on the input
% shift Matrix - dwie redefinicje:
% 1) S := S' (zeby byc "on the same page" z artykulem Pawla)
% 2) squeezed_tilde := -squeezed_tilde:
%    LHS: wektor s jest od biezacego c-sys do innego
%    RHS: wektor s jest od jakiegos c-sys (srodka masy) do biezacego
if length(q) ~= length(input{1})
    error('funkcja getData przyjmuje tylko wektor wsp ZLACZOWYCH');
end

I = eye(3); Om = [0, -1; 1, 0]; 
[fi, omega] = in2outAngles(q, dq, input);

data = struct('M1',[], 'M2',[], 'S12',[], 'S21',[], 'S1c',[],'S2c',[],...
             'Sc1',[],'Sc2',[], 's12',[], 's21',[], 's1C',[],'s2C',[],...
             'sC1',[],'sC2',[],'dS12',[],'dS1c',[],'dS2c',[],...
            'dS21',[],'dSc1',[],'dSc2',[],'H',[],'D',[],'fi',[],...
            'dM1',[],'dM2',[],'dH',[],'dD',[],'om',[]);
N = length(input{1});
for i = 1 : N
    [S12, S21, S1c, Sc1, S2c, Sc2] = deal(I);
    [dS12, dS21, dS1c, dSc1, dS2c, dSc2] = deal(zeros(3));
    M = input{1}(i).M;          sC1_loc = input{1}(i).sC1;
    % note: przyjmujemy s12_loc jako wektor od interfejsu 1 do interfejsu 2
    % (natomiast nie do konca rozwazanego czlunu). Zakladamy, ze do parent
    % body przylaczony jest tylko jeden child body
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
    
    S12(3,1:2) = (Om*s12).';    S21(3,1:2) = (Om*s21).';
    Sc1(3,1:2) = (Om*sC1).';    S1c(3,1:2) = (Om*s1C).';    
    Sc2(3,1:2) = (Om*sC2).';    S2c(3,1:2) = (Om*s2C).';    
    
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
    
    M1 = S1c * M * S1c.';        M2 = S2c * M * S2c.';    
    dM1 = dS1c * M * S1c.' + S1c * M * dS1c.';
    dM2 = dS2c * M * S2c.' + S2c * M * dS2c.';
    
    if strcmp(input{1}(i).type, 'trans')
        u = Rot(fi(i))*[1;0];    v = Om * u;
        H = [u; 0];              D = [ [v;0], [0;0;1] ];
        dH = [Om*u*om; 0];       dD = [[-u*om;0], zeros(3,1)];
    else
        H = input{1}(i).H;       D = input{1}(i).D;
        dH = zeros(3,1);         dD = zeros(3,2);
    end
    
    data(i) = struct('M1',M1,'M2',M2,'S12',S12,'S21',S21,...
               'S1c',S1c,'S2c',S2c,'Sc1',Sc1,'Sc2',Sc2,...
               's12',s12,'s21',s21,'s1C',s1C,'s2C',s2C,...
               'sC1',sC1,'sC2',sC2,...
               'dS12',dS12,'dS1c',dS1c,'dS2c',dS2c,...
               'dS21',dS21,'dSc1',dSc1,'dSc2',dSc2,...
               'H',H,'D',D,'fi',fi(i),'dM1',dM1,'dM2',dM2,...
               'dH',dH,'dD',dD,'om',omega(i));
end
end

