function propagationMap(E, t, lbd, lambda_low, lambda_high, L, spectralScale)
dt = t(2) -t(1);
[slices, nn] = size(E);
Zs = linspace(0,L, slices);
I = abs(E).^2;
if nn > 2^10
    tint = linspace(min(t),max(t),2.^10);
    lbd_int = linspace(lambda_low, lambda_high,2.^10);
    for jk = 1:slices
        spec(jk,:) = interp1(lbd,abs(fftshift(fft(E(jk,:))*dt)/1e-12).^2,lbd_int,'spline');
        Inorm(jk,:) = interp1(t,I(jk,:)./max(I(jk,:)),tint,'spline');
        Snorm(jk,:) = spec(jk,:)./max(spec(jk,:));
    end
    t = tint;
    lbd = lbd_int;
else
    for jk = 1:slices
        spec(jk,:) = (abs(fftshift(fft(E(jk,:))*dt)/1e-12).^2);
        Inorm(jk,:) = (I(jk,:)./max(I(jk,:)));
        Snorm(jk,:) = spec(jk,:)./max(spec(jk,:));
    end
end

figure()
subplot(2,1,1)
pcolor( Zs,t*1e12, Inorm');
xlim([0 L])
ylim([min(t)*1e12 max(t)*1e12])
% view(2)
shading interp
ylabel('Delay (ps)')
colormap jet
if strcmp(spectralScale,'linear')
    subplot(2,1,2)
    pcolor(Zs,lbd*1e9,Snorm');
    ylim([lambda_low*1e9 lambda_high*1e9])
    xlim([0 L])
    ylabel('Wavelength (nm)')
    xlabel('Propagation distance (m)')
    % view(2)
    shading interp
    colormap jet
elseif strcmp(spectralScale,'log')
    subplot(2,1,2)
    pcolor(Zs,lbd*1e9,10*log10(Snorm'+1));
    ylim([lambda_low*1e9 lambda_high*1e9])
    xlim([0 L])
    ylabel('Wavelength (nm)')
    xlabel('Propagation distance (m)')
    % view(2)
    shading interp
    colormap jet
end