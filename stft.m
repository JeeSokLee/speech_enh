function X=stft(x,nfft,nfrm,nshi,wintype)

%STFT: Short-time Fourier transform
%  usage: X=STFT(x,nfft,nfrm,nshi,wintype),

dur  = length(x);
iter = floor((dur-nfrm)/nshi)+1;
w    = window(wintype,nfrm);
X    = zeros(nfft/2+1,iter);

for ind=1:iter
	ibgn     = (ind-1)*nshi;
	q        = x(ibgn+(1:nfrm),1).*w;
	Q        = fft(q,nfft);
	X(:,ind) = Q(1:nfft/2+1,1);
end

end