function G=specsubt(X,N,delta)

%SPECSUBT: spectral (power) subtraction gain estimation
%    usage: G=specsubt(X,N,delta);

G = max(sqrt(max(X-N,0)./X),delta);