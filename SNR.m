function snr=SNR(x,y)
snr = 20 * log10(norm(y(:)) / norm(y(:) - x(:)));
end
