function [Ecomp, Ltot, Dtot] = compressPulse(E, t, f, l0, lc, betas_init, Linit)
format long g

AC_init = autocoTrace(E);
E2 = E;
AC_rad_init = gaussRadius(t,AC_init,'FWHM');
Flag = 0;
Ltot = 0;
while Flag == 0
    E = linearProp(E2,f,l0,lc,betas_init,Linit);
    AC = autocoTrace(E);
    AC_rad = gaussRadius(t, AC, 'FWHM');
    
    if AC_rad < AC_rad_init
        E2 = E;
        AC_rad_init = AC_rad;
        Ltot = Ltot + Linit;
    else
        Linit = Linit/10;
        if Linit < 1e-6;
            Flag = 1;
        end
    end
end

Ecomp = E;
Dtot = Ltot .* betas_init(1);