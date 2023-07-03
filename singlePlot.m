function singlePlot(E, t, lbd, lambda_low, lambda_high, spectralScale)

% This function display the input optical field in both time & wavelength
% domain.
% INPUTS : 
%       E : Complex electric field vector
%       t : Simulation time vector [s]
%       lbd : Simulation wavelength vector [m]
%       lambda_low : Lowest wavelength to display [m]
%       lambda_high : Highest wavelength to display [m]
%       spectralScale : 'linear' for linear scaling & 'log' for dB scaling
%                       in the spectral domain [string]

dt = t(2)-t(1);
res = length(t);
df = 1/res/dt;
f = (-res/2:(res/2)-1).*df;
w = 2*pi.*f;
c = 299792458;
Tspan = max(t);
spec = fftshift(fft(E)*dt)/1e-12;
subplot(2,1,1)
yyaxis left
plot(t*1e12, abs(E).^2)
title(["Pulse width : " (num2str(1e12*2*gaussRadius(t, abs(E).^2,'FWHM'))) "ps. AC width : " num2str(1e12*2*gaussRadius(t, autocoTrace(E),'FWHM')) " ps"])
xlabel('Delay (ps)')
ylabel('Power (W)')
xlim([-Tspan*1e12 Tspan*1e12])
yyaxis right
plot(t*1e12, autocoTrace(E))
ylabel('SHG intensity (u.a)')
if strcmp(spectralScale,'log')
    subplot(2,1,2)
    yyaxis left
    plot(lbd*1e9, 10*log10(abs(spec).^2))
    title(["Spectral width : " num2str(1e9*2*gaussRadius(lbd, abs(spec).^2,'FWHM')) "nm"])
    xlabel('Wavelength (nm)')
    ylabel('Power (dB)')
    xlim([lambda_low*1e9 lambda_high*1e9])
    %yyaxis right
    %plot(lbd*1e9, -angle(E))
    %ylabel('Phase (rad)')
elseif strcmp(spectralScale,'linear')
    subplot(2,1,2)
    %yyaxis left
    plot(lbd*1e9,abs(spec).^2)
    title(["Spectral width : " num2str(1e9*2*gaussRadius(-lbd, abs(spec).^2,'FWHM')) "nm"])
    xlabel('Wavelength (nm)')
    ylabel('Power (W)')
    xlim([lambda_low*1e9 lambda_high*1e9])
    %yyaxis right
    %plot(lbd*1e9, unwrap(-angle(E)))
    %ylabel('Phase (rad)')
else
    fprintf(1,'Unknown scale option\n');
end