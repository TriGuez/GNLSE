function TE = RK4IPGain(E, h, alpha, betas, gamma, fR, hR_w, tau_shock, lbd, wshift, g)
%  Runge-Kutta 4 in interaction picture algorithm for the intelligent adaptative
%  stepsize solver for the generalised nonlinear Schrödinger equation
%  INPUTS : 
%        E : Complex enveloppe of the input optical pulse
%        h : Propagation length [m]
%        alpha : Confinement losses of the fiber [1/m]
%        betas : Taylor coefficients of the propagation constant [s^n./m, n >=2]
%        gamma : Nonlinear coefficient of the fiber [1/W/m]
%        fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
%        in fused silica) []
%        hR_w : Raman scattering response of the fiber in frequency domain
%        tau_shock : Shock time for self-steepening [s]
%        lbd: Wavelength vector of the simulation [m]
%        wshift : Fourier shift of the simulation angular frequency vector 
%        [rad/s] 
%  OUTPUTS : 
%        TE : Complex enveloppe of the output pulse

lbd_l = lbd;
lbd_l(isnan(lbd)) = 0;
alpha = fftshift(silicaLosses(lbd_l)+alpha);
g = fftshift(g);
beta = betas(1)./2*wshift.^2;
for jl = 2:length(betas)
    beta = beta + betas(jl)/factorial(jl+1) * wshift.^(jl+1);
end
opdisphalf = exp((-alpha/2-1i*beta + g./2)*h/2);
TEip = fft(E).*opdisphalf;

if fR > 0
    k1 = opdisphalf.*nonlinearStepFull(E,h,fR,hR_w,gamma,tau_shock,wshift);
    
    Ehalf2 = ifft(TEip+k1/2);
    k2 = nonlinearStepFull(Ehalf2,h,fR,hR_w,gamma,tau_shock,wshift);
    
    Ehalf3 = ifft(TEip+k2/2);
    k3 = nonlinearStepFull(Ehalf3,h,fR,hR_w,gamma,tau_shock,wshift);
    
    Ehalf4 = ifft(opdisphalf.*(TEip+k3));
    k4 = nonlinearStepFull(Ehalf4,h,fR,hR_w,gamma,tau_shock,wshift);
else
    k1 = opdisphalf.*nonlinearStep(E,h,gamma,tau_shock,wshift);
    
    Ehalf2 = ifft(TEip+k1/2);
    k2 = nonlinearStep(Ehalf2,h,gamma,tau_shock,wshift);
    
    Ehalf3 = ifft(TEip+k2/2);
    k3 = nonlinearStep(Ehalf3,h,gamma,tau_shock,wshift);
    
    Ehalf4 = ifft(opdisphalf.*(TEip+k3));
    k4 = nonlinearStep(Ehalf4,h,gamma,tau_shock,wshift);
end
TE = ifft(opdisphalf.*(TEip+k1/6+k2/3+k3/3)+k4/6);