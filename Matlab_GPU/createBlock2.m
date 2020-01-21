function block = createBlock2(fi,p, pnext,input)
M = input.M;                sC2_loc = input.sC2_loc;
sC1_loc = input.sC1_loc;    s12_loc = sC1_loc-sC2_loc;
Om = [0, -1; 1, 0];         H = [0; 0; 1];

s12 = Rot(fi)*s12_loc;      s21 = -Rot(fi)*s12_loc;
s1C = -Rot(fi)*sC1_loc;     s2C = -Rot(fi)*sC2_loc;

S12 = eye(3);               S21 = eye(3);
S1c = eye(3);               S2c = eye(3);
S12(3,1:2) = (Om*s12).';    S21(3,1:2) = (Om*s21).';
S1c(3,1:2) = (Om*s1C).';    S2c(3,1:2) = (Om*s2C).'; 


M1 = S1c * M * S1c.';       M2 = S2c * M * S2c.';

KSI11 = inv(M1);    KSI12 = M1 \ S12;
KSI22 = inv(M2);    KSI21 = M2 \ S21;

KSI10 = M1 \ (+H*p   - S12*H*pnext);
KSI20 = M2 \ (-H*pnext + S21*H*p);

block = [KSI11; KSI12; KSI21; KSI22; KSI10.'; KSI20.'; S12];

function R = Rot(fi)
R = [cos(fi) -sin(fi);
     sin(fi)  cos(fi)];