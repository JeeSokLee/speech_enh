# speech_enhancement
clear;clc;close all

%% settings
wavPath  = '/timit_test_192';
wavFiles = dir([pwd wavPath '/*.wav']);
numFiles = 1;
NumBits = 16;
fs=16000;
nfft = 1024;
nfrm = 512;
nshi = 128;
floor = 0.0001;
snrdb = 10;

numFiles = length(wavFiles);
for iFile=1:numFiles
  iFileName = [pwd, wavPath, ...
                '/' wavFiles(iFile).name(1 : end - 4)];
  clean = audioread([iFileName '.wav']); 
dur = length(clean);

%% gen-noise
noise1 = wavread('white.wav');
noise2 = wavread('babble.wav');
noise3 = wavread('pink.wav');
snr= 10.^(snrdb/20);

[input1,noise1] = gen_noisy(clean,noise1,snrdb);
[input2,noise2] = gen_noisy(clean,noise2,snrdb);
[input3,noise3] = gen_noisy(clean,noise3,snrdb);

%% test
S = stft(clean,nfrm,nfrm,nshi,@hann);
X1 = stft(input1,nfrm,nfrm,nshi,@hann);
ipsd1 = abs(X1) .^ 2;
npsd1 = imcra(ipsd1);

X2 = stft(input2,nfrm,nfrm,nshi,@hann);
ipsd2 = abs(X2) .^ 2;
npsd2 = imcra(ipsd2);

X3 = stft(input3,nfrm,nfrm,nshi,@hann);
ipsd3 = abs(X3) .^ 2;
npsd3 = imcra(ipsd3);

%%
snr_input=SNR(input1,clean);
 
%spectral subtraction
G_ss1 = specsubt(ipsd1,npsd1,floor);
S_ss1 = G_ss1 .* S;
s_ss1 = istft(S_ss1,nfft,nfrm,nshi,@hann);
Y_ss1 = G_ss1 .* X1;
y_ss1 = istft(Y_ss1,nfft,nfrm,nshi,@hann);

G_ss2 = specsubt(ipsd2,npsd2,floor);
S_ss2 = G_ss2 .* S;
s_ss2 = istft(S_ss2,nfft,nfrm,nshi,@hann);
Y_ss2 = G_ss2 .* X2;
y_ss2 = istft(Y_ss2,nfft,nfrm,nshi,@hann);

G_ss3 = specsubt(ipsd3,npsd3,floor);
S_ss3 = G_ss3 .* S;
s_ss3 = istft(S_ss3,nfft,nfrm,nshi,@hann);
Y_ss3 = G_ss3 .* X3;
y_ss3 = istft(Y_ss3,nfft,nfrm,nshi,@hann);

%weiner filter
G_wf1 = wf_enh(ipsd1,npsd1,floor);
S_wf1 = G_wf1 .* S;
s_wf1 = istft(S_wf1,nfft,nfrm,nshi,@hann);
Y_wf1 = G_wf1 .* X1;
y_wf1 = istft(Y_wf1,nfft,nfrm,nshi,@hann);

G_wf2 = wf_enh(ipsd2,npsd2,floor);
S_wf2 = G_wf2 .* S;
s_wf2 = istft(S_wf2,nfft,nfrm,nshi,@hann);
Y_wf2 = G_wf2 .* X2;
y_wf2 = istft(Y_wf2,nfft,nfrm,nshi,@hann);

G_wf3 = wf_enh(ipsd3,npsd3,floor);
S_wf3 = G_wf3 .* S;
s_wf3 = istft(S_wf3,nfft,nfrm,nshi,@hann);
Y_wf3 = G_wf3 .* X3;
y_wf3 = istft(Y_wf3,nfft,nfrm,nshi,@hann);

%maximum likelihood
G_ml1 = mlee(ipsd1,npsd1,floor,snr,0,0.8);
S_ml1 = G_ml1 .* S;
s_ml1 = istft(S_ml1,nfft,nfrm,nshi,@hann);
Y_ml1 = G_ml1 .* X1;
y_ml1 = istft(Y_ml1,nfft,nfrm,nshi,@hann);

G_ml2 = mlee(ipsd2,npsd2,floor,snr,0,0.8);
S_ml2 = G_ml2 .* S;
s_ml2 = istft(S_ml2,nfft,nfrm,nshi,@hann);
Y_ml2 = G_ml2 .* X2;
y_ml2 = istft(Y_ml2,nfft,nfrm,nshi,@hann);

G_ml3 = mlee(ipsd3,npsd3,floor,snr,0,0.8);
S_ml3 = G_ml3 .* S;
s_ml3 = istft(S_ml3,nfft,nfrm,nshi,@hann);
Y_ml3 = G_ml3 .* X3;
y_ml3 = istft(Y_ml3,nfft,nfrm,nshi,@hann);

%imcra
G_omlsa1 = omlsa(ipsd1,npsd1,floor,0.92);
S_omlsa1 = G_omlsa1 .* S;
s_omlsa1 = istft(S_omlsa1,nfrm,nfrm,nshi,@hann);
Y_omlsa1 = G_omlsa1 .* X1;
y_omlsa1 = istft(Y_omlsa1,nfrm,nfrm,nshi,@hann);  

G_omlsa2 = omlsa(ipsd2,npsd2,floor,0.92);
S_omlsa2 = G_omlsa2 .* S;
s_omlsa2 = istft(S_omlsa2,nfrm,nfrm,nshi,@hann);
Y_omlsa2 = G_omlsa2 .* X2;
y_omlsa2 = istft(Y_omlsa2,nfrm,nfrm,nshi,@hann); 

G_omlsa3 = omlsa(ipsd3,npsd3,floor,0.92);
S_omlsa3 = G_omlsa3 .* S;
s_omlsa3 = istft(S_omlsa3,nfrm,nfrm,nshi,@hann);
Y_omlsa3 = G_omlsa3 .* X1;
y_omlsa3 = istft(Y_omlsa3,nfrm,nfrm,nshi,@hann); 

%% SNR
snr_ss1(1,iFile) = SNR(y_ss1,s_ss1);
snr_wf1(1,iFile) = SNR(y_wf1,s_wf1);
snr_ml1(1,iFile) = SNR(y_ml1,s_ml1);
snr_imcra1(1,iFile) = SNR(y_omlsa1,s_omlsa1);

snr_ss2(1,iFile) = SNR(y_ss2,s_ss2);
snr_wf2(1,iFile) = SNR(y_wf2,s_wf2);
snr_ml2(1,iFile) = SNR(y_ml2,s_ml2);
snr_imcra2(1,iFile) = SNR(y_omlsa2,s_omlsa2);

snr_ss3(1,iFile) = SNR(y_ss3,s_ss3);
snr_wf3(1,iFile) = SNR(y_wf3,s_wf3);
snr_ml3(1,iFile) = SNR(y_ml3,s_ml3);
snr_imcra3(1,iFile) = SNR(y_omlsa3,s_omlsa3);

end

%% SNR mean calculate
snr_mean_ss1 = mean(snr_ss1);
snr_mean_wf1 = mean(snr_wf1);
snr_mean_ml1 = mean(snr_ml1);
snr_mean_imcra1 = mean(snr_imcra1);

snr_mean_ss2 = mean(snr_ss2);
snr_mean_wf2 = mean(snr_wf2);
snr_mean_ml2 = mean(snr_ml2);
snr_mean_imcra2 = mean(snr_imcra2);

snr_mean_ss3 = mean(snr_ss3);
snr_mean_wf3 = mean(snr_wf3);
snr_mean_ml3 = mean(snr_ml3);
snr_mean_imcra3 = mean(snr_imcra3);

%% wav write
% wavwrite(clean,16e3,16,'clean');
% wavwrite(input1,16e3,16,'__in_whi');
% wavwrite(input2,16e3,16,'__in_bab');
% wavwrite(input3,16e3,16,'__in_pink');
% 
% wavwrite(y_ss1/2,16e3,16,'_whi_ss');
% wavwrite(y_wf1/2,16e3,16,'_whi_wf');
% wavwrite(y_ml1/2,16e3,16,'_whi_ml');
% wavwrite(y_omlsa1,16e3,16,'_whi_omlsa');
% 
% wavwrite(y_ss2,16e3,16,'_bab_ss');
% wavwrite(y_wf2,16e3,16,'_bab_wf');
% wavwrite(y_ml2,16e3,16,'_bab_ml');
% wavwrite(y_omlsa2,16e3,16,'_bab_omlsa');
% 
% wavwrite(y_ss3,16e3,16,'_pink_ss');
% wavwrite(y_wf3,16e3,16,'_pink_wf');
% wavwrite(y_ml3,16e3,16,'_pink_ml');
% wavwrite(y_omlsa3,16e3,16,'_pink_omlsa');

%% figure plot
figure;
subplot(221);
[B,f,T] = specgram(clean,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Clean speech');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(222);
[B,f,T] = specgram(input1,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Noisy speech (White noise)');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(223);
[B,f,T] = specgram(input2,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Noisy speech (Babble noise)');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(224);
[B,f,T] = specgram(input3,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Noisy speech (Pink noise)');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');


figure;
[B,f,T] = specgram(input1,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B)));axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Noisy speech (White noise)');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
figure;
subplot(221);
[B,f,T] = specgram(y_ss1,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Spectral subtraction');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(222);
[B,f,T] = specgram(y_wf1,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B)));axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Weiner filter');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(223);
[B,f,T] = specgram(y_ml1,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Maximum likelihood');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(224);
[B,f,T] = specgram(y_omlsa1,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B)));axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('IMCRA');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');


figure;
[B,f,T] = specgram(input2,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Noisy speech (Babble noise)');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
figure;
subplot(221);
[B,f,T] = specgram(y_ss2,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Spectral subtraction');
xlabel('Time (sec)');ylabel('Frequency (Hz)');
subplot(222);
[B,f,T] = specgram(y_wf2,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Weiner filter');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(223);
[B,f,T] = specgram(y_ml2,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Maximum likelihood');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(224);
[B,f,T] = specgram(y_omlsa2,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('IMCRA');
xlabel('Time (sec)');ylabel('Frequency (Hz)');


figure;
[B,f,T] = specgram(input3,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Noisy speech (Pink noise)');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
figure;
subplot(221);
[B,f,T] = specgram(y_ss3,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Spectral subtraction');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(222);
[B,f,T] = specgram(y_wf3,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Weiner filter');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(223);
[B,f,T] = specgram(y_ml3,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('Maximum likelihood');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
subplot(224);
[B,f,T] = specgram(y_omlsa3,nfft,fs,hanning(nfrm),nshi);
imagesc(T,f,20*log10(abs(B))); axis xy;
h=colorbar;
set(h, 'ylim', [-80 20]);
title('IMCRA');
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
