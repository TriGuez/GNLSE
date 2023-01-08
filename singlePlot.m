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
Tspan = max(t);
spec = fftshift(fft(E)*dt)/1e-12;
subplot(2,1,1)
plot(t*1e12, abs(E).^2)
xlabel('Delay (ps)')
ylabel('Power (W)')
xlim([-Tspan*1e12 Tspan*1e12])
if strcmp(spectralScale,'log')
    subplot(2,1,2)
    plot(lbd*1e9, 10*log10(abs(spec).^2))
    xlabel('Wavelength (nm)')
    ylabel('Power (dB)')
    xlim([lambda_low*1e9 lambda_high*1e9])
elseif strcmp(spectralScale,'linear')
    subplot(2,1,2)
    plot(lbd*1e9,abs(spec).^2)
    xlabel('Wavelength (nm)')
    ylabel('Power (W)')
    xlim([lambda_low*1e9 lambda_high*1e9])
else
    fprintf(1,'Unknown scale option\n');
end