function [t, dt, f, df, w, lbd, res, Pulse] = reconstructPulse(filename, tspan, llow, lhigh, l0, Ep)

[t, dt, f, df, w, lbd, res] = initGNLSE(tspan, l0, llow);
datas = load(filename);
S = interp1(datas(:,1)*1e-9, datas(:,2),lbd);
S(isnan(S))=0;
S = fftshift(S);
Temp = fftshift(fft(fliplr(S)));
Temp = Temp./sqrt(trapz(t,abs(Temp).^2));
Temp = Temp.*sqrt(Ep);
Pulse = Temp;
singlePlot(Pulse, t, lbd, llow, lhigh, 'linear')