function z=istft(Y,nfft,nfrm,nshi,wintype)

%ISTFT: Inverse short-time Fourier transform
%    usage: z=ISTFT(Y,nfft,nfrm,nshi,wintype);

iter = size(Y,2);
K    = size(Y,1);
dur  = (iter-1)*nshi+nfrm;
z    = zeros(dur,1);
w    = window(wintype,nfrm);

for ind=1:iter
	ibgn = (ind-1)*nshi;
	X    = Y(:,ind);
	X(1) = X(1)/2;
	X(K) = X(K)/2;
	x    = 2*real(ifft(X,nfft));
	z(ibgn+(1:nfrm),1) = z(ibgn+(1:nfrm),1)+x(1:nfrm,1);
end

% for indice=1:iter %Overlapp add technique
%     left_index=((indice-1)*nshi) ;
%     index=left_index+[1:nfrm];
%     temp_ifft=real(ifft(Y(:,indice),nfft));
%     z(index)= z(index)+temp_ifft(index).*w;
% end


end