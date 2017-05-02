function N=imcra(X)

%IMCRA - Improved minima controlled recursive averaging noise estimation
%
%    usage - N = IMCRA(X);


[df,dt] = size(X);

wsub = [0.25;0.5;0.25];
N = zeros(df,dt);

U = 8;
V = 15;
D = U*V;

alpha   = 0.92;
alpha_s = 0.9;
alpha_d = 0.85;
beta    = 1.47;
Bmin    = 1.66;
gamma_0 = 4.6;
gamma_1 = 3;
gamma_2 = gamma_1-1;
zeta_0  = 1.67;

sapmax  = 0.95;

S_hat_sub_min = zeros(df,U);
S_bar_sub_min = zeros(df,U);

%% initialization - first frame
N(:,1) = X(:,1);
N(:,2) = X(:,1);
S_f    = conv(X(:,1),wsub,'same');
S_hat  = S_f;
S_bar  = S_f;

S_hat_sub = S_f;
S_bar_sub = S_f;
for jnd=1:U
	S_hat_sub_min(:,jnd) = S_f;
	S_bar_sub_min(:,jnd) = S_f;
end

S_hat_min = S_f;
S_bar_min = S_f;

Gi         = ones(df,1);
gamma_prev = ones(df,1);
npsd_bias  = N(:,1);
	

%% noise estimation
for tind=2:dt-1 %time index
	ind = rem(tind-1,V)+1;
	S_f = conv(X(:,tind),wsub,'same');
	
	%% SNR related params
	gamma = X(:,tind)./N(:,tind);
	xi    = alpha*(Gi.^2).*gamma_prev+(1-alpha)*max(gamma-1,0);
	nu    = gamma.*xi./(1+xi);
	Gi    = xi./(1+xi).*exp(expint(nu)/2);
	
	%% first smoothing
	S_hat = alpha_s*S_hat+(1-alpha_s)*S_f; % first iteration: gives rough VAD decision
	
	%% first minimum
	if(ind==1)
		S_hat_sub      = S_hat;
	else
		S_hat_sub      = min(S_hat,S_hat_sub);
	end
	S_hat_min          = min(S_hat,S_hat_min);
	S_hat_sub_min(:,1) = S_hat_sub;
	
	%% rough vad
	vad_hat   = zeros(df,1); % indication function
	gamma_hat = X(:,tind)./S_hat_min/Bmin; %posteriori SNR
	zeta_hat  = S_hat./S_hat_min/Bmin; %priori SNR
	vad_hat( gamma_hat<gamma_0 & zeta_hat<zeta_0 ) = 1; %speech is absent
	
	%% second smoothing
	num     = conv(X(:,tind).*vad_hat,wsub,'same');
	den     = conv(vad_hat,wsub,'same');
	S_bar_f = num./den; %smoothing in frequency
	S_bar_f(den==0) = S_bar(den==0);
	S_bar   = alpha_s*S_bar+(1-alpha_s)*S_bar_f; %smoothing in time
  %second iteration: excluded high power speech regions to reduce the variance of minimum
	
	%% second minimum
	if(ind==1)
		S_bar_sub      = S_bar;
	else
		S_bar_sub      = min(S_bar,S_bar_sub);
	end
	S_bar_min          = min(S_bar,S_bar_min);
	S_bar_sub_min(:,1) = S_bar_sub;
	
	%% speech absent/present probability
	gamma_bar = X(:,tind)./S_bar_min/Bmin;
	zeta_bar  = S_hat./S_bar_min/Bmin;
	sap       = (gamma_1-gamma_bar)/gamma_2; %actual speech absence prob
  %gamma_1(threshold) set to satisfy a significant margin so the prob of wrong decision was smaller than a threshold 
	sap(zeta_bar>=zeta_0) = 0;
	sap       = min(max(sap,0),sapmax);
	spp = 1./(1+sap./(1-sap).*(1+xi).*exp(-nu));
	
	%% noise PSD estimation
	alpha_tilde_d = alpha_d+(1-alpha_d)*spp; % time-varying smoothing para.
	npsd_bias     = alpha_tilde_d.*npsd_bias+(1-alpha_tilde_d).*X(:,tind); %update noise spectrum estimate
	N(:,tind+1)   = beta*npsd_bias;
	
	%% finalization
	if(ind==V)
		S_hat_sub_min(:,2:U) = S_hat_sub_min(:,1:U-1);
		S_hat_min            = min(S_hat_sub_min,[],2);
		S_bar_sub_min(:,2:U) = S_bar_sub_min(:,1:U-1);
		S_bar_min            = min(S_bar_sub_min,[],2);
	end
	
	gamma_prev = gamma;
end
