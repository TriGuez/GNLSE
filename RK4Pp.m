function Pzh = RK4Pp(sigma_ap, sigma_ep, h, n2z,Nions, Pp, gamma_P, lbd_p)

k1 =  ((sigma_ep + sigma_ap)*n2z - sigma_ap*Nions)*Pp*gamma_P - (silicaLosses(lbd_p)*Pp);
k2 =  ((sigma_ep + sigma_ap)*n2z - sigma_ap*Nions)*gamma_P*(Pp + (h/2).*k1) - (silicaLosses(lbd_p)*(Pp+(h/2)*k1));
k3 =  ((sigma_ep + sigma_ap)*n2z - sigma_ap*Nions)*gamma_P*(Pp+(h/2)*k2)- (silicaLosses(lbd_p)*(Pp+(h/2)*k2));
k4 =  ((sigma_ep + sigma_ap)*n2z - sigma_ap*Nions)*gamma_P*(Pp+h*k3)- (silicaLosses(lbd_p)*(Pp+(h)*k3));

Pzh = Pp + h*((k1/6)+(k2/3)+(k3/3)+(k4/6));