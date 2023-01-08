function hR_w = ramanResponseBW(fR, wshift)
%  This function computes the Raman scattering response of an optical fiber in
%  the frequency domain according to Blow, K. J., & Wood, D. (1989). 
%  Theoretical description of transient stimulated Raman scattering in optical 
%  fibers. IEEE Journal of Quantum Electronics, 25(12), 2665–2673.
%  INPUTS : 
%        fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
%        in fused silica) []
%        wshift : Fourier shift of the simulation angular frequency vector 
%        [rad/s] 
%  OUTPUTS:
%        hR_w : Raman scattering response in the frequency domain

t1 = 12.2e-15;
t2 = 32e-15;
if fR > 0
    hR_w = conj((t1^2 + t2^2)./(t2^2-t1^2.*(1i+t2*wshift).^2));
else
    hR_w = zeros(1,length(wshift));
end