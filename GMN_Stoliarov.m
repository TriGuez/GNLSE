clear
close all
format long g
c = 299792458;
Tspan = 10e-12;
lambda_low = 980e-9;
lambda_high = 1200e-9;
l0 = 1100e-9;
[t, dt, f, df, w, lbd, res] = initGNLSE(Tspan, l0, lambda_low);
f0 = c./l0;
% Pulse definition
P = 40;
tFWHM_p = 0.5e-12;
t0_p = tFWHM_p/2/sqrt(log(2));
C2 = 0;
t_shift =8e-12;
l1 = 1035e-9;
f1 = c./l1;
frep = 100e6;
E = gaussianPulse(P, C2, t0_p, f1,t_shift, t, f, f0);
singlePlot(E, t, lbd, lambda_low, lambda_high, 'linear')

% Fiber definition
betas = [1.8583e-26, 3.8775e-41, -4.7658e-56, 1.6251e-70];
gamma = 2.4705e-3;
L = 5;
slices = 100;
h = L/slices;
r_core = (11.5/2)*1e-6;
r_clad = (125/2)*1e-6;
fR = 0.18;
Gamma_P = (pi*r_core.^2)./(pi*r_clad.^2);
[sigma_a, sigma_e] = VLMAcrossSections(lbd);
N_ions = 7.45113843884923e25;
Pp = 5;
lbd_p = 976e-9;
[sigma_ap, sigma_ep] = VLMAcrossSections(lbd_p);
tol = 1e-5;
[Eout,~,~,pout] = propagationFibreGain(E, h, h/2, l0, 1040e-9, tol, t, f, lbd, 0,...
    betas, gamma,fR, frep,Pp, lbd_p,sigma_a,sigma_e,sigma_ap, sigma_ep,...
    N_ions, r_core,Gamma_P,'LMA-10-125');
Esave = Eout';
for jk = 2:slices
    [Eout,~,~,pout] = propagationFibreGain(Eout, h, h/2, l0, 1040e-9, tol, t, f, lbd, 0,...
    betas, gamma,fR, frep,pout, lbd_p,sigma_a,sigma_e,sigma_ap, sigma_ep,...
    N_ions, r_core,Gamma_P,'LMA-10-125');
    Esave = [Esave Eout'];
end
Esave = Esave';
propagationMap(Esave,t,lbd,lambda_low,lambda_high,L,'log')
bcomp = GratingCompressor(1200,45,l0);
figure()
[Ecomp, Ltot] = compressPulse(Eout,t,f,l0,l0,bcomp,1e-4,lambda_low,1300e-9,lbd);