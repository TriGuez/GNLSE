clear
close all

c = 299792458;
Tspan = 300e-12;


% Signal band
lc_sig = 810e-9;
llow_sig = 800e-9;
f0_sig = c./lc_sig;
w0_sig = 2*pi*f0_sig;
[t, dt, f_sig, df, w_sig, lbd_sig, res] = initGNLSE(Tspan, lc_sig, llow_sig);
Esig = gaussianPulse(10000,1.3e-23,105e-15,f0_sig,0,t,f_sig,f0_sig);

% Pump band
lc_p = 1040e-9;
f_p = (-res/2:res/2-1).*df;
f0_p = c./lc_p;
w0_p = 2*pi*f0_p;
lbd_p = c./(f_p+f0_p);
Epump = gaussianPulse(10000,1.3e-23,105e-15,f0_p,0,t,f_p,f0_p);


% Idler band
w0_i = w0_p - (w0_sig - w0_p);
f0_i = w0_i./(2*pi);
f_i = (-res/2:res/2-1).*df;
lbd_i = c./(f_i+f0_i);
Eidl = gaussianPulse(10000,1.3e-23,105e-15,f0_i,0,t,f_i,f0_i);

lbd_full = [lbd_sig; lbd_p; lbd_i];
U = [Esig; Epump; Eidl];

plotBNLSE(t, lbd_full, U)