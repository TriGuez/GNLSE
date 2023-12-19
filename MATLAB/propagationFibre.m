function [Eout, temp_rad, spec_rad] = propagationFibre(E, L, h, l0, lc, tol ,t, f, lbd, alpha, betas, gamma, fR, FiberName)

% This function computes the propagation over a specific distance of an
% optical fiber.
% INPUTS : 
%       E : Complex enveloppe of the input electrical field
%       L : Length of the optical fiber [m]
%       h : Propagation initial stepsize [m]
%       l0 : Central vacuum wavelength of the simulation [m]
%       lc : Central vacuum wavelength of the input optical field [m]
%       tol : Relative tolerance of the simulation []
%       t : Simulation time vector [s]
%       f : Simulation frequency vector [Hz]
%       lbd : Simulation wavelength vector [m]
%       alpha : Confinement losses of the fiber [1/m]
%       betas : Taylor coefficients of the propagation constant @ lc [s^n/m
%               n >=2]
%       gamma : Nonlinear coefficient of the fiber [1/W/m]
%       fR : Fractionnal Raman response of the propagation medium (0.18 in
%       fused silica). Set to zero to exclude Raman scattering []
%       FiberName : String to display while computing
% OUTPUTS :
%       Eout : Complex enveloppe of the output optical field
%       temp_rad : pulse duration evolution during the propagation [s]
%       spec_rad : Pulse spectral width evolution during the propagation

c = 299792458;
f0 = c/l0;
fc = c/lc;
f = (f+f0) - fc;
wshift = fftshift(2*pi*f);
hR_w = ramanResponseBW(fR, wshift);
tau_shock = 1/(2*pi*f0);
[Eout, temp_rad, spec_rad] = adaptativeSolver(E, L, h, alpha, betas,...
            gamma, fR, hR_w, tau_shock, t, lbd, wshift, tol, FiberName);