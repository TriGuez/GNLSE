clear
close all

c = 299792458;
Tspan = 600e-12;
lambda_low = 700e-9;
lambda_high = 1700e-9;
l0 = 1040e-9;
f0 = c./l0;
[t, dt, f, df, w, lbd, res, lambda_low, lambda_high] = initGNLSE(Tspan, l0, lambda_low, lambda_high);
wshift = fftshift(w);

J_p = 3e-6;
tFWHM_p = 100e-15;
t0_p = tFWHM_p/2/sqrt(log(2));
C2_p = 1.3e-23;
l1 = l0;
l2 = 802e-9;
f2 = c./l2;
f_p = c/l1;
tshift = 0;
w0 = 2*pi*f_p;
P_p = 0.94*(J_p/tFWHM_p);
Epump = gaussianPulse(P_p,C2_p,t0_p,f_p, tshift, t, f, f0);
singlePlot(Epump, t, lbd, lambda_low, lambda_high, 'log')
spec = fftshift(fft(Epump)*dt)/1e-12;
phi_f = angle(spec);
phi_t = angle(Epump);

fR = 0.18;
% 
gamma_SUP5 = 0.012;
alpha = 0;
betas_SUP5 = eval_beta('neff_SUP5.txt',l0,2,12);

L_SUP5 = 8e-2;
h = L_SUP5/200;
tol = 1e-5;
Esig = rectPulse(1e-3,Tspan-1e-12, f2, t, f, f0);
     E = Eout+Epump;
     [Eamp, ~, ~] = propagationFibre(E, L_SUP5, h, l0, l1, tol, t, f, lbd,...
                 alpha, betas_SUP5, gamma_SUP5, fR, 'SUP5');
     singlePlot(Eout, t, lbd, lambda_low, lambda_high, 'linear');
