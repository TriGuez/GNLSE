function [Ecomp, Ltot] = compressPulse(E, t, f, l0, lc, betas, Linit, lambda_low, lambda_high, lbd)
format long g
E = centerPulse(E,t);
AC1 = autocoTrace(E);
rad1 = gaussRadius(t, AC1, '1/e2');

E2 = linearProp(E, f, l0, lc, betas, Linit);
AC2 = autocoTrace(E2);
rad2 = gaussRadius(t, AC2, '1/e2');
Ltot = Linit;

if rad2 < rad1
    while rad2<rad1
        Ltot = Ltot+Linit;
        E = E2;
        AC1 = autocoTrace(E);
        rad1 = gaussRadius(t, AC1, '1/e2');
        E2 = linearProp(E, f, l0, lc, betas, Linit);
        drawnow
        AC2 = autocoTrace(E2);
        rad2 = gaussRadius(t, AC2, '1/e2');
        singlePlot(E2, t, lbd, lambda_low, lambda_high, 'linear')
        subplot(2,1,1)
        % xlim([-(rad2+rad2*3) (rad2+rad2*3)]*1e12)
    end
    Ecomp = E;
    singlePlot(E, t, lbd, lambda_low, lambda_high, 'linear')
    subplot(2,1,1)
    % xlim([-(rad2+rad2*3) (rad2+rad2*3)]*1e12)
else
    Ltot = 0;
    Ecomp = 0;
end
 