function k = nonlinearStepFull(E, h, fR, hR_w, gamma, tau_shock, wshift)
%  Nonlinear quarter-step for the RK4IP algorithm including Raman scattering
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
op1 = fR.*E.*ifft(hR_w.*fft(abs(E).^2));
op2 = fft(E.*((1-fR)*abs(E).^2)+op1);
k = -h.*1i.*gamma.*(1+wshift.*tau_shock).*op2;