% Exemple de la propagation et du décalage Raman d'un soliton d'ordre 2
% dans une fibre avec et sans animation

clear
close all
set(0,'DefaultFigureWindowStyle','docked')

% Définition des paramètres de simus & des axes utiles
c = 299792458;
Tspan = 4e-12;
lambda_low = 900e-9;
lambda_high = 1300e-9;
l0 = 1040e-9;
f0 = c./l0;
[t, dt, f, df, w, lbd, res, lambda_low, lambda_high] = initGNLSE(Tspan, l0, lambda_low, lambda_high);
wshift = fftshift(w);
tol = 1e-7;

% Définition des paramètres de la fibre
alpha = 0;
betas = [-3.53571099317077e-26, 3.68095336905838e-41,2.0172694409917e-55...
         ,-1.31835263681886e-69];
gamma = 0.0742400655658342;
L = 1;
fR = 0.18;
% Définition de l'impulsion initiale

N = 2;
tFWHM_p = 100e-15;
t0_p = tFWHM_p/2/sqrt(log(2));
l1 = l0;
f_p = c/l1;
tshift = -3e-12;
P_p = (N^2*abs(betas(1)))./(t0_p^2*gamma);
Epump = sechPulse(P_p,0,t0_p,f_p,tshift,t,f,f0);

% Paramètres pour animation + sauvegarde dans une matrice
slices = 500;
dL = L/slices;
h = dL/1000;
Esave = zeros(slices+1,length(t));
Eout = Epump;
Esave(1,:) = Epump;
for jk = 2:slices
    Eout = propagationFibre(Eout, dL, h, l0, l0, tol, t, f, lbd,...
             alpha, betas, gamma, fR, 'Propagation example');
    singlePlot(Eout, t, lbd, lambda_low, lambda_high, 'linear')
    title(['L = ' num2str(jk*dL) ' m'])
    Esave(jk,:) = Eout;
    drawnow
end

% Simulation sans animation
% Eout = propagationFibre(Epump,L,h,l0,l1,tol,t,f,lbd,alpha,betas,gamma,...
%     fR, 'Raman self frequency shift');
% singlePlot(Eout, t, lbd, lambda_low, lambda_high, 'linear')
propagationMap(Esave(1:end,:),t,lbd,lambda_low,lambda_high,L,'linear');
