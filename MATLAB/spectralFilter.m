function Eout = spectralFilter(E, l_c, lFWHM, t, lbd, filterType)
%  Computes a spectral filtering of the input optical field
%  INPUTS : 
%        E : Complex enveloppe of the input optical field
%        l_c : Central wavelength of the filter [m]
%        lFWHM : filter's full width at half maximum [m]
%        t : Time vector [s]
%        lbd : Wavelength vector [m]
%        filterType : Shape of the filter. Default is 'gaussian'
%  OUTPUTS : 
%        Eout : Complex enveloppe of the filtered optical field

dt = t(2)-t(1);
if strcmp(filterType, 'gaussian')
    eL = lFWHM/2/sqrt(log(2));
    filtre = exp(-0.5*((lbd-l_c)./eL).^2);
elseif strcmp(filterType, 'square')
    filtre = zeros(1,length(lbd));
    filtre((lbd > (l_c - lFWHM)) & (lbd < (l_c + lFWHM))) = 1;
elseif strcmp(filterType, 'LPF')
    filtre = zeros(1,length(lbd));
    filtre(lbd < l_c) = 1;
elseif strcmp(filterType, 'HPF')
    filtre = zeros(1, length(lbd));
    filtre(lbd > l_c) = 1;
else
    fprintf(1,'Unknown filter shape\n');
end
filtre(isnan(filtre)) = 0;
spec = (fftshift(fft(E)*dt)./1e-12).*filtre;
Eout =ifft(fftshift(spec).*1e-12)./dt;
