function isEqual = meq(v1, v2, eps)
% przyrównuje do siebie dwie wielkoœci w oparciu o dok³adnoœæ eps
if nargin == 2
    eps = 1e-10; end

if all(abs(v1 - v2) <= eps)
    isEqual = true;
else
    isEqual = false;
end