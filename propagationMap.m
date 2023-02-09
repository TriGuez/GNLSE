function propagationMap(E, t, lbd, lambda_low, lambda_high, L)
N = 10;
dt = t(2) -t(1);
[slices, ~] = size(E);
Zs = linspace(0,L, slices);
I = abs(E).^2;
tint = linspace(min(t),max(t),2.^N);
lbd_int = linspace(lambda_low, lambda_high,2.^N);
for jk = 1:slices
    spec(jk,:) = interp1(lbd,abs(fftshift(fft(E(jk,:))*dt)/1e-12).^2,lbd_int,'spline');
    Inorm(jk,:) = interp1(t,I(jk,:)./max(I(jk,:)),tint,'spline');
    Snorm(jk,:) = spec(jk,:)./max(spec(jk,:));
end
t = tint;
lbd = lbd_int;
figure()
subplot(2,1,1)
pcolor( Zs,t*1e12, Inorm');
xlim([0 L])
ylim([min(t)*1e12 max(t)*1e12])
% view(2)
shading interp
ylabel('Delay (ps)')
colormap jet

subplot(2,1,2)
pcolor(Zs,lbd*1e9,Snorm');
ylim([lambda_low*1e9 lambda_high*1e9])
xlim([0 L])
ylabel('Wavelength (nm)')
xlabel('Propagation distance (m)')
% view(2)
shading interp
colormap jet
