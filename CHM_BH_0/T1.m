function Y=T1(X)
% 20200519
% shortcut for pTrans(...)
s = size(X,1);
d = sqrt(s); % square (x) only!
Y=Gamma(X,1,[d d]);

