function k = nonlinearStep(E, h, gamma, tau_shock, wshift)
%  Nonlinear quarter-step for the RK4IP algorithm without Raman scattering
%  INPUTS : 
%        E : Complex enveloppe of the input optical pulse
%        h : Propagation length [m]
%        fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
%        in fused silica) []
%        hR_w : Raman scattering response of the fiber in frequency domain
%        gamma : Nonlinear coefficient of the fiber [W.?¹.m?¹]
%        tau_shock : Shock time for self-steepening [s]
%        wshift : Fourier shift of the simulation angular frequency vector 
%        [rad.s?¹] 
%  OUTPUTS : 
%        k : Nonlinear quarter-step

op1 = fft(E.*abs(E).^2);
k = -h.*1i.*gamma.*(1+wshift.*tau_shock).*op1;