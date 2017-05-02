function G=omlsa(X,N,delta,alpha)

%OMLSA: Optimally-modified log-spectral amplitude gain estimator
%
%    usage: G=omlsa(X,N,delta,alpha);

wlocal    = 1;
wglobal   = 15;

hlocal    = window(@hann,2*wlocal+1);
hlocal    = hlocal/sum(hlocal);
hglobal   = window(@hann,2*wglobal+1);
hglobal   = hglobal/sum(hglobal);

zetamin   = 10^(-10/10);
zetamax   = 10^(-5/10);
zetaratio = log(zetamax/zetamin);

zetapmin  = 10^(0/10);
zetapmax  = 10^(10/10);

sapmax    = 0.95;
beta      = 0.7;

[df,dt] = size(X);

G = zeros(df,dt);
for ind=1:dt
	npsd  = N(:,ind);
	xpsd  = X(:,ind);
	
	%% snr related params estimation
	gamma = xpsd./npsd;
	if(ind==1)
		xi = alpha+(1-alpha)*max(gamma-1,0);
	else
		xi = alpha*xi_old+(1-alpha)*max(gamma-1,0);
	end
	nu = xi./(1+xi).*gamma;
	
	%% sap/spp related params estimation
	if ind==1
		zeta      = xi;
		zetaframe = mean(zeta);
		zetapeak  = zetaframe;
		sap       = sapmax*ones(df,1);
	else
		zeta       = beta*zetaprev+(1-beta)*xi_prev;
		zetalocal  = conv(zeta,hlocal,'same');
		zetaglobal = conv(zeta,hglobal,'same');
		zetaframe  = mean(zeta);

		mu      = max(min(log(zetaframe/zetapeak/zetamin)/zetaratio,1),0);
		Plocal  = max(min(log(zetalocal/zetamin)/zetaratio,1),0);
		Pglobal = max(min(log(zetaglobal/zetamin)/zetaratio,1),0);
		if(zetaframe>zetamin)
			if(zetaframe>zetaframeprev)
				Pframe = 1;
				zetapeak = min(max(zetaframe,zetapmin),zetapmax);
			else
				Pframe = mu;
			end
		else
			Pframe = 0;
		end
		
		sap = min(1-Pframe*Plocal.*Pglobal,sapmax);
	end
	spp = 1./(1+sap./(1-sap).*(1+xi).*exp(-nu));
	
	%% gain estimation
	Gi       = xi./(1+xi).*exp(expint(nu)/2);
	G(:,ind) = (Gi.^spp).*(delta.^(1-spp));
	
	%% finalization
	xi_old        = (Gi.^2).*gamma;
	xi_prev       = xi;
	zetaframeprev = zetaframe;
	zetaprev      = zeta;
end


















