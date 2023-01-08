function AC = autocoTrace(E)

if ~isreal(E)
    E = abs(E).^2;
end
AC = fftshift(ifft(fft(E).*conj(fft(E))));
AC = AC./max(AC);



