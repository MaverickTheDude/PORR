function data = getS(q, input)
% returns auxaliary matrices and vector for current state based on the input
% shift Matrix - dwie redefinicje:
% 1) S := S' (zeby byc "on the same page" z artykulem Pawla)
% 2) squeezed_tilde := -squeezed_tilde:
%    LHS: wektor s jest od biezacego c-sys do innego
%    RHS: wektor s jest od jakiegos c-sys (srodka masy) do biezacego
if length(q) ~= length(input{1})
    error('funkcja getS przyjmuje tylko wektor wsp ZLACZOWYCH');
end
I = eye(3); Om = [0, -1; 1, 0]; 
fi = in2outAngles(q, input);

data = struct('M1',[], 'M2',[], 'S12',[], 'S21',[], 'S1c',[],'S2c',[],...
             'Sc1',[],'Sc2',[], 's12',[], 's21',[], 's1C',[],'s2C',[],...
             'sC1',[],'sC2',[], 'H',[], 'D',[], 'fi',[]);
N = length(input{1});
for i = 1 : N
    [S12, S21, S1c, Sc1, S2c, Sc2] = deal(I);
    M = input{1}(i).M;          sC1_loc = input{1}(i).sC1;
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
    S1c(3,1:2) = (Om*s1C).';    Sc1(3,1:2) = (Om*sC1).';
    S2c(3,1:2) = (Om*s2C).';    Sc2(3,1:2) = (Om*sC2).';
    
    M1 = S1c * M * S1c.';       M2 = S2c * M * S2c.';
    
    if strcmp(input{1}(i).type, 'trans')
        u = Rot(fi(i))*[1;0];    v = Om * u;
        H = [u; 0];              D = [ [v;0], [0;0;1] ];
    else
        H = input{1}(i).H;       D = input{1}(i).D;
    end
    
    data(i) = struct('M1',M1,'M2',M2,'S12',S12,'S21',S21,...
               'S1c',S1c,'S2c',S2c,'Sc1',Sc1,'Sc2',Sc2,...
               's12',s12,'s21',s21,'s1C',s1C,'s2C',s2C,...
               'sC1',sC1,'sC2',sC2,'H',H,'D',D,'fi',fi(i));
end
end