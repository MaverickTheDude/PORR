function KSI = ksiCoefs(q, p, input)
% i in name refers to "index"
bodyNo = length(q);
data = getS(q, input);
KSI = struct('i11',[],'i12',[],'i21',[],'i22',[],'i10',[],'i20',[]);

for i = 1:bodyNo-1
    M1 = data(i).M1;    M2 = data(i).M2;
    S12 = data(i).S12;  S21 = data(i).S21;
    KSI11 = inv(M1);    KSI12 = M1 \ S12;
    KSI22 = inv(M2);    KSI21 = M2 \ S21;

    KSI10 = M1 \ (+data(i).H*p(i)   - S12*data(i+1).H*p(i+1));
    KSI20 = M2 \ (-data(i+1).H*p(i+1) + S21*data(i).H*p(i));
    KSI(i) = struct('i11',KSI11,'i12',KSI12,'i21',KSI21,'i22',KSI22,...
                 'i10',KSI10,'i20',KSI20); 
end
% Terminal body
M1 = data(bodyNo).M1;    M2 = data(bodyNo).M2;
S12 = data(bodyNo).S12;  S21 = data(bodyNo).S21;
KSI11 = inv(M1);    KSI12 = M1 \ S12;
KSI22 = inv(M2);    KSI21 = M2 \ S21;

KSI10 = M1 \ (data(bodyNo).H*p(bodyNo)) ;
KSI20 = M2 \ (S21*data(bodyNo).H*p(bodyNo));
KSI(bodyNo) = struct('i11',KSI11,'i12',KSI12,'i21',KSI21,'i22',KSI22,...
    'i10',KSI10,'i20',KSI20);