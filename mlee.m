function G=mlee(X,N,delta,snr,suflag,q)

% MLEE : maximum likelihood envelope estimation
% g=MLEE(x,n,delta,snr,suflag,q);

xi = 10^(snr/20);
Gi = max((1+sqrt(max(X-N,0)./X))/2,delta);
if(suflag)
	p = min(2*sqrt(xi.*X./N),400);
	P = exp(-xi)*besseli(0,p);
	Q = (q*P)./(1+(q*P));
    %Q = (P)./(1+(P));
else
	Q = 1;
end
G = Gi.*Q;
	
end