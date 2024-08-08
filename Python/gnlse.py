import numpy as np
import pyfftw
import multiprocessing
import tqdm
import matplotlib.pyplot as plt
from misc_math import *
from math import factorial


def initGNLSE(Tspan, l0, lambda_low, lambda_high):
# This function creates the usefull vectors & values 
# INPUTS : 
#       Tspan : half width of the temporal window [s]
#       l0 : Central wavelength of the simulation [m]
#       lambda_low : lowest wavelength of the simulation [m]
# OUTPUTS : 
#       t : Time vector from -Tspan to Tspan [s]
#       dt : Time resolution [s]
#       f : Frequency vector [Hz]
#       df : Frequency resolution [Hz]
#       lbd : Wavelength vector [m]
#       res : Size of the simulation []
    c = 299792458 # vacuum speed of light [m.s**-1]
    f0 = c/l0
    Fspan = c/lambda_low - f0
    nu = 1

    while (4*Tspan*Fspan > 2**nu) :
        nu = nu +1
    res = 2**nu
    dt = 2*Tspan/(res-1)
    t = np.arange(-res/2,res/2)*dt
    df = 1/res/dt
    f = np.arange(-res/2,res/2)*df
    w = 2*np.pi*f
    lbd = c/(f+f0)
    lbd[lbd<0] = float('nan')

    idx_range = np.argwhere((lbd>lambda_low) & (lbd < lambda_high))
    if ((np.isnan(lambda_low)) & np.isnan(lambda_high)) :
        idx_range = np.argwhere(np.isnan(lbd))
        idx_range = np.asarray([np.max(idx_range)+1, np.size(lbd)-1])
        lambda_low = np.min(lbd)
        lambda_low = np.max(lbd)
    else :
        if (lambda_low < lbd[idx_range[np.size(idx_range)-1]]) :
            lambda_low = np.round(lbd[idx_range[np.size(idx_range)-1]]*1e9)*1e-9
        if (lambda_high > lbd[idx_range[0]]) :
            lambda_high = np.round(lbd[idx_range[0]]*1e9)*1e-9

    return t, dt, f, df, w, lbd, res, lambda_low, lambda_high

def spectralFilter(E, l_c, lFWHM, t, lbd, filterType='gaussian') : 
# Computes a spectral filtering of the input optical field
# INPUTS : 
#       E : Complex enveloppe of the input optical field
#       l_c : Central wavelength of the filter [m]
#       lFWHM : filter's full width at half maximum [m]
#       t : Time vector [s]
#       lbd : Wavelength vector [m]
#       filterType : Shape of the filter. Default is 'gaussian'
# OUTPUTS : 
#       Eout : Complex enveloppe of the filtered optical field
    dt = t[1]-t[0]
    if filterType == 'gaussian' :
        eL = lFWHM/2/np.sqrt(np.log(2))
        filtre = np.exp(-0.5*((lbd-l_c)/eL)**2)
    if filterType == 'rect' : 
        filtre = np.zeros(len(lbd))
        filtre[np.argwhere((lbd > (l_c-lFWHM)) & (lbd < (l_c + lFWHM)))] = 1
    filtre[np.isnan(filtre)] = 0
    spec_out = (np.fft.fftshift(np.fft.fft(E)*dt)/1e-12)*filtre
    Eout = np.fft.ifft(np.fft.fftshift(spec_out)*1e-12)/dt
    return Eout, spec_out


def gaussianPulse(P, C2, t0, f1, tshift, t, f, f0):
# This function computes the complex enveloppe of a gaussian optical pulse
# centered at frequency f1 and optical noise of one photon per spectral node
# INPUTS : 
#       P : Peak power of the fourier transform limited pulse [W]
#       C2 : 2nd order spectral phase of the pulse [sÂ²]
#       t0 : 1/e radius of the fourier transform limited pulse [s]
#       f1 : Central frequency of the pulse [Hz]
#       tshift : Time offset of the pulse [s]
#       t : Time vector of the simulation [s]
#       f : Frequency vector of the pulse [Hz]
#       f0 : Central frequency of the simulation [Hz]
# OUTPUTS :
#       E : Complex enveloppe of the optical pulse 
    dt = t[1]-t[0]
    df = f[1] - f[0]
    res = np.size(t)
    h = 6.62607004e-34
    wshift = np.fft.fftshift(2*np.pi*f)

    noise = np.sqrt(h*np.fft.fftshift(f+f0)/df*(res-1)/res)*np.exp(-1j*np.random.rand(res,)*2*np.pi)
    noise[np.isnan(noise)]=0+0*1j
    spectral_phase = C2/2*wshift**2
    E = np.sqrt(P)*np.exp(-0.5*(((t-tshift)/t0)**2))*np.exp(-1j*2*np.pi*(f0-f1)*t)
    E = np.fft.ifft(np.fft.fft(E)*dt*np.exp(-1j*spectral_phase)+(noise))/dt
    return E

def rectPulse(P, tFWHM, f1, t, f, f0):
# This function computes the complex enveloppe of a rectangular optical pulse
# centered at frequency f1 and optical noise of one photon per spectral node
# INPUTS : 
#       P : Peak power of the pulse [W]
#       tFWHM : Half temporal width of the pulse [s]
#       f1 : Central frequency of the pulse [Hz]
#       t : Time vector of the simulation [s]
#       f : Frequency vector of the pulse [Hz]
#       f0 : Central frequency of the simulation [Hz]
# OUTPUTS :
#       E : Complex enveloppe of the optical pulse 
    dt = t[1]-t[0]
    df = f[1] - f[0]
    res = len(t)
    h = 6.62607004e-34
    noise = np.sqrt(h*np.fft.fftshift(f+f0)/df*(res-1)/res)*np.exp(-1j*np.random.rand(res,)*2*np.pi)
    noise[np.isnan(noise)]=0+0*1j
    spectral_phase = np.zeros((np.size(t)))
    E = np.sqrt(P)*((t>-tFWHM) & (t<tFWHM))*np.exp(-1j*2*np.pi*(f0-f1)*t)
    E = np.fft.ifft(np.fft.fft(E)*dt*np.exp(-1j*spectral_phase)+noise)/dt
    return E

def sechPulse(P, C2, t0, f1, tshift, t, f, f0) :
# This function computes the complex enveloppe of a gaussian optical pulse
# centered at frequency f1 and optical noise of one photon per spectral node
# INPUTS : 
#       P : Peak power of the fourier transform limited pulse [W]
#       C2 : 2nd order spectral phase of the pulse [sÂ²]
#       t0 : 1/e radius of the fourier transform limited pulse [s]
#       f1 : Central frequency of the pulse [Hz]
#       tshift : Time offset of the pulse [s]
#       t : Time vector of the simulation [s]
#       f : Frequency vector of the pulse [Hz]
#       f0 : Central frequency of the simulation [Hz]
# OUTPUTS :
#       E : Complex enveloppe of the optical pulse 
    dt = t[1]-t[0]
    df = f[1] - f[0]
    res = np.size(t)
    h = 6.62607004e-34
    wshift = np.fft.fftshift(2*np.pi*f)

    noise = np.sqrt(h*np.fft.fftshift(f+f0)/df*(res-1)/res)*np.exp(-1j*np.random.rand(res,)*2*np.pi)
    noise[np.isnan(noise)]=0+0*1j
    spectral_phase = C2/2*wshift**2
    E = np.sqrt(P)*sech((t-tshift)/t0)*np.exp(-1j*2*np.pi*(f0-f1)*t)
    E = np.fft.ifft(np.fft.fft(E)*dt*np.exp(-1j*spectral_phase)+(noise))/dt
    return E

def autocoTrace(E) :
    x = pyfftw.empty_aligned(np.size(E), dtype = 'complex128')
    X = pyfftw.empty_aligned(np.size(E), dtype = 'complex128')
    fft = pyfftw.FFTW(x,X,threads = multiprocessing.cpu_count())
    ifft = pyfftw.FFTW(X,x,direction = 'FFTW_BACKWARD', threads = multiprocessing.cpu_count())
    x[:] = E*E
    convo = fft()
    X[:] = convo
    return ifft()

def ramanResponseBW(fR, wshift) :
# This function computes the Raman scattering response of an optical fiber in
# the frequency domain according to Blow, K. J., & Wood, D. (1989). 
# Theoretical description of transient stimulated Raman scattering in optical 
# fibers. IEEE Journal of Quantum Electronics, 25(12), 2665â€2673.
# INPUTS : 
#       fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
#       in fused silica) []
#       wshift : Fourier shift of the simulation angular frequency vector 
#       [rad.s?Â¹] 
# OUTPUTS:
#       hR_w : Raman scattering response in the frequency domain
    t1 = 12.2e-15
    t2 = 32e-15
    if fR > 0 :
        hR_w = np.conj((t1**2+t2**2)/(t2**2-t1**2*(1j+t2*wshift)**2))
    else :
        hR_w = np.zeros(np.size(wshift))
    hR_w[np.isnan(hR_w)]=0
    return hR_w


def silicaLosses(lbd) :
# This function computes the silica optical losses vs wavelength, according to:
# Sorensen, S. T. " Deep-blue supercontinuum light sources based on tapered
# photonic crystal fibers", PhD thesis, DTU Fotonik ,2013.
# INPUTS :
#       lbd : Wavelength vector [m]
# OUTPUTS : 
#       losses : Losses vector [m?Â¹]
    auv = 0.001*np.exp((4.67e-6)/lbd)/1000/4.343
    air = 6e11*np.exp(-(47.8e-6)/lbd)/1000/4.343
    asc = ((1.3/((lbd*1e6)**4))+1)/1000/4.343
    aoh = (7/(1+((lbd-1.380e-6)/16e-9)**2))/1000/4.343
    losses = auv + air + asc + aoh
    for j in range(0,np.size(lbd)):
        if lbd[j]<300e-9 :
            losses[j] = 0
        if lbd[j]>2500e-9 :
            losses[j] = 0
    return losses


def nonlinearStepFull(E, h, fR, hR_w, gamma, tau_shock, wshift) : 
# Nonlinear quarter-step for the RK4IP algorithm including Raman scattering
# INPUTS : 
#       E : Complex enveloppe of the input optical pulse
#       h : Propagation length [m]
#       fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
#       in fused silica) []
#       hR_w : Raman scattering response of the fiber in frequency domain
#       gamma : Nonlinear coefficient of the fiber [W.?Â¹.m?Â¹]
#       tau_shock : Shock time for self-steepening [s]
#       wshift : Fourier shift of the simulation angular frequency vector 
#       [rad.s?Â¹] 
# OUTPUTS : 
#       k : Nonlinear quarter-step
    x = pyfftw.empty_aligned(np.size(wshift), dtype = 'complex128')
    X = pyfftw.empty_aligned(np.size(wshift), dtype = 'complex128')
    fft = pyfftw.FFTW(x,X, threads = multiprocessing.cpu_count())
    ifft = pyfftw.FFTW(X,x,direction = 'FFTW_BACKWARD', threads = multiprocessing.cpu_count())
    x[:] = np.abs(E)**2
    op1 = fft()
    X[:] = hR_w*op1
    op1 = fR*E*ifft()
    x[:] = E*((1-fR)*np.abs(E)**2)+op1
    op2 = fft()
    k = -h*1j*gamma*(1+wshift*tau_shock)*op2
    return k


def silicaLosses(lbd) :
# This function computes the silica optical losses vs wavelength, according to:
# Sorensen, S. T. " Deep-blue supercontinuum light sources based on tapered
# photonic crystal fibers", PhD thesis, DTU Fotonik ,2013.
# INPUTS :
#       lbd : Wavelength vector [m]
# OUTPUTS : 
#       losses : Losses vector [m?Â¹]

    auv = 0.001*np.exp((4.67e-6)/lbd)/1000/4.343
    air = 6e11*np.exp(-(47.8e-6)/lbd)/1000/4.343
    asc = ((1.3/((lbd*1e6)**4))+1)/1000/4.343
    aoh = (7/(1+((lbd-1.380e-6)/16e-9)**2))/1000/4.343
    losses = auv + air + asc + aoh
    if np.size(lbd) > 1 :
        for j in range(0,np.size(lbd)):
            if lbd[j]<300e-9 :
                losses[j] = 0
            if lbd[j]>2500e-9 :
                losses[j] = 0
    return losses


def nonlinearStepFull(E, h, fR, hR_w, gamma, tau_shock, wshift) : 
# Nonlinear quarter-step for the RK4IP algorithm including Raman scattering
# INPUTS : 
#       E : Complex enveloppe of the input optical pulse
#       h : Propagation length [m]
#       fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
#       in fused silica) []
#       hR_w : Raman scattering response of the fiber in frequency domain
#       gamma : Nonlinear coefficient of the fiber [W.?Â¹.m?Â¹]
#       tau_shock : Shock time for self-steepening [s]
#       wshift : Fourier shift of the simulation angular frequency vector 
#       [rad.s?Â¹] 
# OUTPUTS : 
#       k : Nonlinear quarter-step
    x = pyfftw.empty_aligned(np.size(wshift), dtype = 'complex128')
    X = pyfftw.empty_aligned(np.size(wshift), dtype = 'complex128')
    fft = pyfftw.FFTW(x,X, threads = multiprocessing.cpu_count())
    ifft = pyfftw.FFTW(X,x,direction = 'FFTW_BACKWARD', threads = multiprocessing.cpu_count())
    x[:] = np.abs(E)**2
    op1 = fft()
    X[:] = hR_w*op1
    op1 = fR*E*ifft()
    x[:] = E*((1-fR)*np.abs(E)**2)+op1
    op2 = fft()
    k = -h*1j*gamma*(1+wshift*tau_shock)*op2
    return k

def nonlinearStep(E, h, gamma, tau_shock, wshift) :
# Nonlinear quarter-step for the RK4IP algorithm without Raman scattering
# INPUTS : 
#       E : Complex enveloppe of the input optical pulse
#       h : Propagation length [m]
#       fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
#       in fused silica) []
#       hR_w : Raman scattering response of the fiber in frequency domain
#       gamma : Nonlinear coefficient of the fiber [W.?Â¹.m?Â¹]
#       tau_shock : Shock time for self-steepening [s]
#       wshift : Fourier shift of the simulation angular frequency vector 
#       [rad.s?Â¹] 
# OUTPUTS : 
#       k : Nonlinear quarter-step
    x = pyfftw.empty_aligned(np.size(wshift), dtype = 'complex128')
    X = pyfftw.empty_aligned(np.size(wshift), dtype = 'complex128')
    fft = pyfftw.FFTW(x,X, threads = multiprocessing.cpu_count())

    x[:] = E*np.abs(E)**2
    op1 = fft()
    k = -h*1j*gamma*(1+wshift*tau_shock)*op1

    return k


def RK4ip(E, h, alpha, betas, gamma, fR, hR_w, tau_shock, lbd, wshift) :
# Runge-Kutta 4 in interaction picture algorithm for the intelligent adaptative
# stepsize solver for the generalised nonlinear SchrÃ¶dinger equation
# INPUTS : 
#       E : Complex enveloppe of the input optical pulse
#       h : Propagation length [m]
#       alpha : Confinement losses of the fiber [m?Â¹]
#       betas : Taylor coefficients of the propagation constant [s^n.m?Â¹]
#       gamma : Nonlinear coefficient of the fiber [W.?Â¹.m?Â¹]
#       fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
#       in fused silica) []
#       hR_w : Raman scattering response of the fiber in frequency domain
#       tau_shock : Shock time for self-steepening [s]
#       lbd: Wavelength vector of the simulation [m]
#       wshift : Fourier shift of the simulation angular frequency vector 
#       [rad.s?Â¹] 
# OUTPUTS : 
#       TE : Complex enveloppe of the output pulse
    y = pyfftw.empty_aligned(np.size(wshift), dtype = 'complex128')
    Y = pyfftw.empty_aligned(np.size(wshift), dtype = 'complex128')
    fft = pyfftw.FFTW(y,Y, threads = multiprocessing.cpu_count())
    ifft = pyfftw.FFTW(Y,y,direction = 'FFTW_BACKWARD', threads = multiprocessing.cpu_count())

    lbd_l = lbd
    lbd_l[np.isnan(lbd_l)]=0
    alpha = np.fft.fftshift(silicaLosses(lbd_l)+alpha)

    beta = betas[0]/2*wshift**2
    for jl in range(1,np.size(betas)) :
        beta = beta + betas[jl]/factorial(jl+2)*wshift**(jl+2)
    
    opdisphalf = np.exp((-alpha/2-1j*beta)*h/2)
    y[:] = E
    TEip = fft()*opdisphalf

    if fR > 0 :
        k1 = opdisphalf*nonlinearStepFull(E,h,fR,hR_w,gamma,tau_shock,wshift)

        Y[:] = TEip+k1/2
        Ehalf2 = ifft()
        k2 = nonlinearStepFull(Ehalf2,h,fR,hR_w,gamma,tau_shock,wshift)

        Y[:] = TEip + k2/2
        Ehalf3 = ifft()
        k3 = nonlinearStepFull(Ehalf3,h,fR,hR_w,gamma,tau_shock,wshift)

        Y[:] = opdisphalf*(TEip + k3)
        Ehalf4 = ifft()
        k4 = nonlinearStepFull(Ehalf4,h,fR,hR_w,gamma,tau_shock,wshift)
    else :
        k1 = opdisphalf*nonlinearStep(E, h, gamma, tau_shock, wshift)

        Y[:] = TEip+k1/2
        Ehalf2 = ifft()
        k2 = nonlinearStep(Ehalf2, h, gamma, tau_shock, wshift)

        Y[:] = TEip+k2/2
        Ehalf3 = ifft()
        k3 = nonlinearStep(Ehalf3, h, gamma, tau_shock, wshift)

        Y[:] = opdisphalf*(TEip+k3)
        Ehalf4 = ifft()
        k4 = nonlinearStep(Ehalf4, h, gamma, tau_shock, wshift)
    
    Y[:] = opdisphalf*(TEip+k1/6+k2/3+k3/3) + k4/6
    return ifft()


def adaptiveSolver(E, L , h, alpha, betas, gamma, fR, hR_w, tau_shock, t, lbd, wshift, tol):
# This function computes the complex enveloppe of an optical pulse during its
# propagation in a given optical fiber, by solving the following generalised 
# nonlinear SchrÃ¶dinger equation for the complex field A : 
# dz(A) = -(alpha(w)/2)A + sum_{k>=2}(i**k+1/k!)beta_k*d(t**k)A) + i*gamma*...
#   ... (1+i*tau_shock*dt)*(A*int_{-inf}^{inf}(R(T')*abs(A(z,T-T'))Â²dT'))
# The equation is solved using a Runge-Kutta 4 in interaction picture (RK4IP)
# (Balac, S. & al. "The interaction picture method for solving the generalized
# nonlinear SchrÃ¶dinger equation in optics", ESAIM: M2AN, vol.50, nÂ°4, p.945-964)
# The solver is an intelligent adaptative stepsize solver, from Nguyen, D. T.
# "Modeling and Design Photonics by Examples using Matlab", IOP publishing, 
# 2021. doi : 10.1088/978-0-7503-2272-0
# INPUTS : 
#       E : Complex enveloppe of the input optical pulse
#       L : Length of the optical fiber [m]
#       h : Initial propagation stepsize [m]
#       alpha : Confinement losses of the fiber [m?Â¹]    n = len(y)
#       betas : Taylor coefficients of the propagation constant [s^n.m?Â¹]
#       gamma : Nonlinear coefficient of the fiber [W.?Â¹.m?Â¹]
#       fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
#       in fused silica) []
#       hR_w : Raman scattering response of the fiber in frequency domain
#       tau_shock : Shock time for self-steepening [s]
#       lbd: Wavelength vector of the simulation [m]
#       wshift : Fourier shift of the simulation angular frequency vector 
#       [rad.s?Â¹] 
#       tol : Maximum tolerance for the adaptative stepsize
# OUTPUTS :
#       Eout : Complex enveloppe of the optical pulse at the output of the fiber
    progress_bar = tqdm.tqdm(total=L,unit='m')
    Z_prop = 0
    while Z_prop < L:
        progress_bar.n = Z_prop
        Z_prop = Z_prop + h
        Uf = RK4ip(E, h, alpha, betas, gamma, fR, hR_w, tau_shock, lbd, wshift)
        Uc = RK4ip(RK4ip(E,h/2,alpha,betas,gamma,fR,hR_w,tau_shock,lbd,wshift),h/2,alpha,betas,gamma,fR,hR_w,tau_shock,lbd,wshift)
        error = np.sqrt(np.sum(np.abs(Uf-Uc)**2))/np.sqrt(np.sum(np.abs(Uf)**2))
        factor = tol/error

        if error > 2*tol :
            Z_prop = Z_prop-h
        else : 
            E = (16/15)*Uf-(1/15)*Uc
        
        h = h*factor**(1/5)

        if Z_prop + h > L :
            h = L - Z_prop
        progress_bar.update(0)
    progress_bar.close()
    return E

def propagationFibre(E, L, h, l0, lc, tol, t, f, lbd, alpha, betas, gamma, fR) :
    c = 299792458
    f0 = c/l0
    fc = c/lc
    f = (f+f0) - fc
    wshift = np.fft.fftshift(2*np.pi*f)
    hR_w = ramanResponseBW(fR,wshift)
    tau_shock = 1/(2*np.pi*fc)
    return adaptiveSolver(E, L, h, alpha, betas, gamma, fR, hR_w, tau_shock, t, lbd, wshift, tol)


def singlePlot(E, t, lbd, lambda_low, lambda_high, spectralScale = 'log') :
    dt = t[1]-t[0]
    Tspan = np.max(t)
    spec = np.fft.fftshift(np.fft.fft(E)*dt)/1e-12
    fig, axs = plt.subplots(2)
    axs[0].plot(t*1e12, 1e-3*np.abs(E)**2)
    axs[0].set_xlim((-Tspan*1e12, Tspan*1e12))
    axs[0].set_xlabel('Delay (ps)')
    axs[0].set_ylabel('Power (kW)')
    if spectralScale == 'log' :
        axs[1].plot(lbd*1e9,10*np.log10(np.abs(spec)**2))
        axs[1].set_ylabel('Power (dB)')
    if spectralScale == 'linear' :
        axs[1].plot(lbd*1e9, np.abs(spec)**2)
        axs[1].set_ylabel('Power (W)')
    axs[1].set_xlim((lambda_low*1e9, lambda_high*1e9))
    axs[1].set_xlabel('Wavelength (nm)')

    plt.show()

def n2Yb(wshift, rho_w, sigma_a, sigma_e, sigma_ap, sigma_ep, lbd_p, Pp, Gamma_P, frep, rcore, Nions) :
    tau = 0.88e-3
    hbar = (6.62607015e-34)/(2*np.pi)
    wp = 2*np.pi*(299792458/lbd_p)
    num = (((Gamma_P*sigma_ap)/(hbar*wp))*Pp) + frep*np.trapz(((sigma_a*rho_w)/(hbar*wshift)), wshift)
    denom = (((Gamma_P*(sigma_ap+sigma_ep))/(hbar*wp))*Pp) + frep*np.trapz((((sigma_a+sigma_e)*rho_w)/(hbar*wshift)), wshift) + ((np.pi*rcore**2)/tau)
    return Nions*(num/denom)

def RK4Pp(sigma_ap, sigma_ep, h, n2z,Nions, Pp, gamma_P, lbd_p) :
    k1 =  ((sigma_ep + sigma_ap)*n2z - sigma_ap*Nions)*Pp*gamma_P - (silicaLosses(lbd_p)*Pp)
    k2 =  ((sigma_ep + sigma_ap)*n2z - sigma_ap*Nions)*gamma_P*(Pp + (h/2)*k1) - (silicaLosses(lbd_p)*(Pp+(h/2)*k1))
    k3 =  ((sigma_ep + sigma_ap)*n2z - sigma_ap*Nions)*gamma_P*(Pp+(h/2)*k2)- (silicaLosses(lbd_p)*(Pp+(h/2)*k2))
    k4 =  ((sigma_ep + sigma_ap)*n2z - sigma_ap*Nions)*gamma_P*(Pp+h*k3)- (silicaLosses(lbd_p)*(Pp+(h)*k3))

    return Pp + h*((k1/6)+(k2/3)+(k3/3)+(k4/6))

def adaptiveSolverGain(E, L , h, alpha, betas, gamma, fR, hR_w, tau_shock, t, lbd, wshift, tol,
                        frep, Pp, lbd_p, sigma_a,sigma_e, sigma_ap, sigma_ep, N_ions, r_core,Gamma_P,ww0):
    progress_bar = tqdm.tqdm(total=L,unit='m')
    Z_prop = 0
    dt = t[1]-t[0]
    Pump_out = 0
    S = np.fft.fftshift(np.fft.fft(E)*dt)/np.sqrt(2*np.pi)
    rho_w = np.abs(S)**2
    while Z_prop < L:
        progress_bar.n = Z_prop
        Z_prop = Z_prop + h
        n_up = n2Yb(ww0, rho_w, sigma_a, sigma_e, sigma_ap, sigma_ep, lbd_p, Pp, Gamma_P, frep, r_core, N_ions)
        Pp = RK4Pp(sigma_ap, sigma_ep, h, n_up, N_ions, Pp, Gamma_P, lbd_p)
        S = np.fft.fftshift(np.fft.fft(E) * dt) / np.sqrt(2 * np.pi)
        rho_w = np.abs(S) ** 2
        g = (sigma_a + sigma_e)*n_up - sigma_a*N_ions
        Uf = RK4IPGain(E, h, alpha, betas, gamma, fR, hR_w, tau_shock, lbd, wshift, g)
        Uc = RK4IPGain(RK4IPGain(E,h/2,alpha,betas,gamma,fR,hR_w,tau_shock,lbd,wshift, g),h/2,alpha,betas,gamma,fR,hR_w,tau_shock,lbd,wshift,g)
        error = np.sqrt(np.sum(np.abs(Uf-Uc)**2))/np.sqrt(np.sum(np.abs(Uf)**2))
        factor = tol/error

        if error > 2*tol :
            Z_prop = Z_prop-h
        else : 
            E = (16/15)*Uf-(1/15)*Uc
        
        h = h*factor**(1/5)

        if Z_prop + h > L :
            h = L - Z_prop
        progress_bar.update(0)
    progress_bar.close()
    return E, Pump_out

def RK4IPGain(E, h, alpha, betas, gamma, fR, hR_w, tau_shock, lbd, wshift,g) :
    y = pyfftw.empty_aligned(np.size(wshift), dtype='complex128')
    Y = pyfftw.empty_aligned(np.size(wshift), dtype='complex128')
    fft = pyfftw.FFTW(y, Y, threads=multiprocessing.cpu_count())
    ifft = pyfftw.FFTW(Y, y, direction='FFTW_BACKWARD', threads=multiprocessing.cpu_count())

    lbd_l = lbd
    lbd_l[np.isnan(lbd_l)] = 0
    alpha = np.fft.fftshift(silicaLosses(lbd_l) + alpha)

    beta = betas[0] / 2 * wshift ** 2
    for jl in range(1, np.size(betas)):
        beta = beta + betas[jl] / factorial(jl + 2) * wshift ** (jl + 2)
    g = np.fft.fftshift(g)
    opdisphalf = np.exp((-alpha / 2 - 1j * beta + g/2) * h / 2)
    y[:] = E
    TEip = fft() * opdisphalf

    if fR > 0:
        k1 = opdisphalf * nonlinearStepFull(E, h, fR, hR_w, gamma, tau_shock, wshift)

        Y[:] = TEip + k1 / 2
        Ehalf2 = ifft()
        k2 = nonlinearStepFull(Ehalf2, h, fR, hR_w, gamma, tau_shock, wshift)

        Y[:] = TEip + k2 / 2
        Ehalf3 = ifft()
        k3 = nonlinearStepFull(Ehalf3, h, fR, hR_w, gamma, tau_shock, wshift)

        Y[:] = opdisphalf * (TEip + k3)
        Ehalf4 = ifft()
        k4 = nonlinearStepFull(Ehalf4, h, fR, hR_w, gamma, tau_shock, wshift)
    else:
        k1 = opdisphalf * nonlinearStep(E, h, gamma, tau_shock, wshift)

        Y[:] = TEip + k1 / 2
        Ehalf2 = ifft()
        k2 = nonlinearStep(Ehalf2, h, gamma, tau_shock, wshift)

        Y[:] = TEip + k2 / 2
        Ehalf3 = ifft()
        k3 = nonlinearStep(Ehalf3, h, gamma, tau_shock, wshift)

        Y[:] = opdisphalf * (TEip + k3)
        Ehalf4 = ifft()
        k4 = nonlinearStep(Ehalf4, h, gamma, tau_shock, wshift)

    Y[:] = opdisphalf * (TEip + k1 / 6 + k2 / 3 + k3 / 3) + k4 / 6
    return ifft()
def dummyCrossYb(lbd) :
    Ai_a = [299.92*1e-27, 2314.71*1e-27, 75.54*1e-27, 211.61*1e-27, 496.62*1e-27, 5.12*1e-27]
    Bi_a = [973.12*1e-9, 977.05*1e-9, 1014.99*1e-9, 908.36*1e-9, 917.03*1e-9, 1059.54*1e-9]
    Ci_a = [32.84*1e-9, 7.24*1e-9, 37.69*1e-9, 92.3*1e-9, 48.67*1e-9, 40.2*1e-9]

    Ai_e = [268.98*1e-27, 633.48*1e-27, 104.74*1e-27, 323.84*1e-27, 2299.03*1e-27]
    Bi_e = [981.68*1e-9, 1025.96*1e-9, 975.32*1e-9, 1072.08*1e-9, 977.53*1e-9]
    Ci_e = [27.16*1e-9, 42.01*1e-9, 86.16*1e-9, 44.96*1e-9, 7.25*1e-9]

    sigma_a = np.zeros(np.size(lbd))
    sigma_e = np.zeros(np.size(lbd))

    for jk in range(np.size(Ai_a)) :
        sigma_a = sigma_a+(Ai_a[jk]*np.exp(-(((lbd-Bi_a[jk]))/Ci_a[jk])**2))
    for jk in range(np.size(Ai_e)) :
        sigma_e = sigma_e + (Ai_e[jk]*np.exp(-(((lbd-Bi_e[jk]))/Ci_e[jk])**2))

    return sigma_a, sigma_e

def propagationFibreGain(E, L, h, l0, lc, tol, t, f, lbd, alpha, betas, gamma, fR, frep, Pp, lbd_p, sigma_a, sigma_e,
                         sigma_ap, sigma_ep, N_ions, r_core, Gamma_P) :

    c = 299792458
    f0 = c/l0
    fc = c/lc
    f = (f+f0) - fc
    ww0 = (2*np.pi*f) + (2*np.pi*fc)
    wshift = np.fft.fftshift(2*np.pi*f)
    hR_w = ramanResponseBW(fR, ww0)
    tau_shock = 1/(2*np.pi*f0)
    Eout, Pump_out = adaptiveSolverGain(E, L , h, alpha, betas, gamma, fR, hR_w, tau_shock, t, lbd, wshift, tol,
                        frep, Pp, lbd_p, sigma_a,sigma_e, sigma_ap, sigma_ep, N_ions, r_core,Gamma_P,ww0)
    return Eout, Pump_out