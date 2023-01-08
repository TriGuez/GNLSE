function [t, dt, f, df, w, lbd, res] = initGNLSE(Tspan, l0, lambda_low)

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