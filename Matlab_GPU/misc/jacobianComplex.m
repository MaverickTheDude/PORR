function Fun_q = jacobianComplex(fun, q0, input)
%verifyDerivatives_Complex uses complex step differential scheme to
%numerically calculate jacobian of the vector function fun. 

h = 1e-8;
n = length(q0);
Nf = length(fun(q0, input));
Fun_q = zeros(Nf,n);
% Complex differences
for i = 1:n
    delta = zeros(n,1); delta(i) = 1i*h;
    funPerturbed = fun(q0 + delta, input);
    Fun_q(:,i) = imag(funPerturbed) / h;  
end

% Real differences
%{
for i = 1:n
    delta = zeros(n,1); delta(i) = h;
    funForward = fun(q0 + delta, input);
    funRev =  fun(q0 - delta, input);
    Fun_q(:,i) = (funForward - funRev) / (2*h);  
end
%}