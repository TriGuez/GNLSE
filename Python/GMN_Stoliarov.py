from gnlse import *

c = 299792458
Tspan = 15e-12
lambda_low = 900e-9
lambda_high = 1300e-9
l0 = 1100e-9

[t, dt, f, df, w, lbd, res, lambda_low, lambda_high] = initGNLSE(Tspan, l0, lambda_low, lambda_high)
f0 = c/l0

P = 40
tFWHM = 1e-12
t0 = tFWHM/2/np.sqrt(np.log(2))
C2 = 0
tshift = 8e-12
l1 = 1040e-9
f1 = c/l1
frep = 10e6
E = gaussianPulse(P, C2, t0, f1, tshift, t, f, f0)
singlePlot(E, t, lbd, lambda_low, lambda_high)

alpha = 0
betas = [1.8583e-26, 3.8775e-41, -4.7658e-56, 1.6251e-70]
gamma = 2.4705e-3
L = 5
h = L/100
r_core = 0.5*(11.5e-6)
r_clad = 0.5*(125e-6)
fR = 0.18
Gamma_P = (np.pi*r_core**2)/(np.pi*r_clad**2)
sigma_a, sigma_e = dummyCrossYb(lbd)
N_ions = 7.45113843884923e25
P_p = 1
lbd_p = 976e-9
sigma_ap, sigma_ep = dummyCrossYb(lbd_p)
tol = 1e-7
Eout, Pump_out = propagationFibreGain(E, L, h, l0, l1, tol , t, f, lbd, alpha, betas, gamma, fR, frep, P_p, lbd_p,
                                      sigma_a, sigma_e, sigma_ap, sigma_ep, N_ions, r_core, Gamma_P)
singlePlot(Eout, t, lbd, lambda_low, lambda_high, 'linear')