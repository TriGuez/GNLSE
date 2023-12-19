function [t, dt, f, df, w, lbd, res, lambda_low, lambda_high] = initGNLSE(Tspan, l0, lambda_low, lambda_high)

% This function creates the usefull vectors & values
% INPUTS : 
%       Tspan : half width of the temporal window [s]
%       l0 : Central vacuum wavelength of the simulation [m]
%       lambda_low : Lowest vacuum wavelength of the simulation [m]
% OUTPUTS : 
%       t : Time vector from -Tspan to Tspan [s]
%       dt : time stepsize [s]
%       f : Relative frequency vector [Hz]
%       df : frequency stepsize [Hz]
%       w : Angular frequency vector [rad/s]
%       lbd : Wavelength vector [m]
%       res : Number of points in the simulation []

c = 299792458;
f0 = c/l0;
Fspan = c/lambda_low - f0;
nu = 1;
while 4*Tspan*Fspan > 2^nu
    nu = nu + 1;
end
res = 2^nu;
dt = 2*Tspan/(res-1);
t = (-res/2:(res/2)-1)*dt;
df = 1/res/dt;
f = (-res/2:(res/2)-1)*df;
w = 2*pi*f;
lbd = c./(f+f0);
lbd(lbd<0) = nan;

idx_range = find((lbd>lambda_low & lbd < lambda_high));
if ( isnan(lambda_low)&&isnan(lambda_high) )
    idx_range = find(isnan(lbd));
    idx_range = [max(idx_range)+1 length(lbd)-1];
    lambda_low = min(lbd);
    lambda_high = max(lbd);
else
    if lambda_low < lbd(idx_range(end))
        lambda_low = round(lbd(idx_range(end))*1e9).*1e-9;
    end
    if lambda_high > lbd(idx_range(1))
        lambda_high = round(lbd(idx_range(1))*1e9).*1e-9;
    end
end