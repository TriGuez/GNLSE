clear
close all
format long g
c = 299792458;
Tspan = 10e-12;
lambda_low = 980e-9;
lambda_high = 1300e-9;
l0 = 1100e-9;
[t, dt, f, df, w, lbd, res] = initGNLSE(Tspan, l0, lambda_low);
f0 = c./l0;
% Pulse definition
P = 1000;
tFWHM_p = 80e-15;
t0_p = tFWHM_p/2/sqrt(log(2));
% t0_p = 200e-15;
C2 = 0;
t_shift = 5e-12;
l1 = 1035e-9;
f1 = c./l1;
frep = 80e6;
E = gaussianPulse(P, C2, t0_p, f1,t_shift, t, f, f0);
singlePlot(E, t, lbd, lambda_low, lambda_high, 'linear')

% Fiber definition
betas = [1.3e-26, 3.9e-41, -4.8e-56, 1.6e-70];
gamma = 1e-3;
L = 5;
h = L/2000;
r_core = (11.5/2)*1e-6;
r_clad = (125/2)*1e-6;
fR = 0.2;
Gamma_P = (pi*r_core.^2)./(pi*r_clad.^2);
[sigma_a, sigma_e] = VLMAcrossSections(lbd);
N_ions = 2.9318e26;
Pp = 2.5;
lbd_p = 976e-9;
[sigma_ap, sigma_ep] = VLMAcrossSections(lbd_p);

Eout = propagationFibreGain(E, L, h, l0, l1, 1e-5, t, f, lbd, 0, betas, gamma,fR, frep,Pp, lbd_p,sigma_a,sigma_e,sigma_ap, sigma_ep, N_ions, r_core,Gamma_P,'test');