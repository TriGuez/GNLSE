from gnlse import *
import numpy as np

c = 299792458
Tspan = 5e-12
lambda_low = 900e-9
lambda_high = 1300e-9
l0 = 1040e-9
f0 = c/l0
[t, dt, f, df, w, lbd, res] = initGNLSE(Tspan, l0, lambda_low)

wshift = np.fft.fftshift(w)
tol = 1e-7

alpha = 0
betas = [-3.53571099317077e-26, 3.68095336905838e-41,2.0172694409917e-55,
         -1.31835263681886e-69]

gamma = 0.0742400655658342
L = 1
fR = 0.18

N=2
tFWHM_p = 100e-15
t0_p = tFWHM_p/2/sqrt(np.log(2))
l1 = l0
f_p = c/l1
tshift = -3e-12
P_p = (N**2*abs(betas[0]))/(t0_p**2*gamma)
Epump = sechPulse(P_p, 0, t0_p, f_p, tshift, t, f, f0)

Eout = propagationFibre(Epump, L, L/200, l0, l1, tol, t, f, lbd, alpha, betas, gamma, fR)

singlePlot(Eout, t, lbd, lambda_low, lambda_high, 'linear')