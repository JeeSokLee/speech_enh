function [Y,N]=gen_noisy(X,N,SNR)

%% Calculate length
lenX = length(X);
lenN = length(N);
len = min(lenX, lenN);

%% Noise signal
if lenX <= lenN
  N = N(1:len, :);
else
  N = [N ; N(lenX-len, :)];
end

%% Power
Ps = mean(diag(X' * X));
Pn = N' * N;
N  = N * sqrt(Ps / Pn);

%% Output
Y = X + (10 ^ (-SNR / 20)) * N;

end
