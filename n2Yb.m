function n_up = n2Yb(wshift ,rho_w,sigma_a,sigma_e, sigma_ap, sigma_ep, lbd_p, Pp, gamma_P, frep, rcore, N_ions)

tau = 88e-3; % 88 ms for Yb
hbar = (6.62607015e-34)/(2*pi); % normalized Planck constant
wp = 2*pi*(299792458/lbd_p); % Pump angular frequency


num = (((gamma_P.*sigma_ap)/(hbar*wp))*Pp) + frep*trapz(wshift,((sigma_a.*rho_w)./(hbar.*wshift)));
denom = (((gamma_P*(sigma_ap+sigma_ep))/(hbar*wp))*Pp) + frep*trapz(wshift,(((sigma_a+sigma_e).*rho_w)./(hbar.*wshift))) + ((pi*rcore^2)./tau);

n_up = N_ions .* (num./denom);