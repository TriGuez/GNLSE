function plotBNLSE(t, lbd, U)

dt = t(2)-t(1);
subplot(2,1,1)
plot(t*1e12, abs(U(1,:)).^2,t*1e12, abs(U(2,:)).^2,t*1e12, abs(U(3,:)).^2)
xlabel('Delay (ps)')
ylabel('Power (W)')
legend('Signal', 'Pump', 'Idler')

subplot(2,3,4)
plot(lbd(1,:)*1e9,abs(fftshift(fft(U(1,:)*dt)./1e-12)).^2)
xlabel('Wavelength (nm)')
ylabel('Power (W)')
legend('Signal')

subplot(2,3,5)
plot(lbd(2,:)*1e9,abs(fftshift(fft(U(2,:)*dt)./1e-12)).^2)
xlabel('Wavelength (nm)')
ylabel('Power (W)')
legend('Pump')

subplot(2,3,6)
plot(lbd(3,:)*1e9,abs(fftshift(fft(U(3,:)*dt)./1e-12)).^2)
xlabel('Wavelength (nm)')
ylabel('Power (W)')
legend('Idler')