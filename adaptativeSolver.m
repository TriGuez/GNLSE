function [Eout, temp_rad, spec_rad] = adaptativeSolver(E, L , h, alpha, betas, gamma, fR,...
            hR_w, tau_shock,t, lbd, wshift, tol, FiberName)
%  This function computes the complex enveloppe of an optical pulse during its
%  propagation in a given optical fiber, by solving the following generalised 
%  nonlinear Schrödinger equation for the complex field A : 
%  dz(A) = -(alpha(w)/2)A + sum_{k>=2}(i**k+1/k!)beta_k*d(t**k)A) + i*gamma*...
%    ... (1+i*tau_shock*dt)*(A*int_{-inf}^{inf}(R(T')*abs(A(z,T-T'))²dT'))
%  The equation is solved using a Runge-Kutta 4 in interaction picture (RK4IP)
%  (Balac, S. & al. "The interaction picture method for solving the generalized
%  nonlinear Schrödinger equation in optics", ESAIM: M2AN, vol.50, n°4, p.945-964)
%  The solver is an intelligent adaptative stepsize solver, from Nguyen, D. T.
%  "Modeling and Design Photonics by Examples using Matlab", IOP publishing, 
%  2021. doi : 10.1088/978-0-7503-2272-0
%  INPUTS : 
%        E : Complex enveloppe of the input optical pulse
%        L : Length of the optical fiber [m]
%        h : Initial propagation stepsize [m]
%        alpha : Confinement losses of the fiber [m?¹]    
%        betas : Taylor coefficients of the propagation constant [s^n.m?¹]
%        gamma : Nonlinear coefficient of the fiber [W.?¹.m?¹]
%        fR : Fractionnal Raman response of the propagation medium (typ. 0.18 
%        in fused silica) []
%        hR_w : Raman scattering response of the fiber in frequency domain
%        tau_shock : Shock time for self-steepening [s]
%        lbd: Wavelength vector of the simulation [m]
%        wshift : Fourier shift of the simulation angular frequency vector 
%        [rad.s?¹] 
%        tol : Maximum tolerance for the adaptative stepsize
%  OUTPUTS :
%        Eout : Complex enveloppe of the optical pulse at the output of the fiber
str = [FiberName ' : '];
textprogressbar(str);
Z_prop = 0;
dt = t(2)-t(1);
temp_rad = [];
spec_rad = [];
epulse = trapz(t, abs(E).^2);
while Z_prop < L
    Z_prop = Z_prop + h;
    Uf = RK4IP(E, h, alpha, betas, gamma, fR, hR_w, tau_shock, lbd, wshift);
    Uc = RK4IP(RK4IP(E, h/2, alpha, betas, gamma, fR, hR_w, tau_shock, lbd,...
            wshift), h/2, alpha, betas, gamma, fR, hR_w, tau_shock, lbd, wshift);
    error = sqrt(sum(abs(Uf-Uc).^2))./sqrt(sum(abs(Uf).^2));
    factor = tol/error;
    
    if error > 2*tol
        Z_prop = Z_prop - h;
    else
        E = (16/15).*Uf - (1/15).*Uc;
        temp_rad = [temp_rad gaussRadius(t,abs(E).^2,'1/e')];
        Espec = fftshift(fft(E)*dt)/1e-12;
        spec_rad = [spec_rad gaussRadius(lbd,abs(Espec).^2,'1/e')];
        textprogressbar((Z_prop/L)*100);
    end
    h = h*factor^(1/5);
    
end
eout = trapz(t,abs(E).^2);
var_phot = (1-(eout/epulse))*100;
textprogressbar(' Done!')
fprintf(1, ['var phot : ' num2str(var_phot) '%% \n']);
Eout = E;