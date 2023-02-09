clear
close all

c = 299792458;
Tspan = 600e-12;
lambda_low = 750e-9;
lambda_high = 1600e-9;
l0 = 1040e-9;
f0 = c./l0;
[t, dt, f, df, w, lbd, res] = initGNLSE(Tspan, l0, lambda_low);
wshift = fftshift(w);

J_p = 1e-6;
tFWHM_p = 105e-15;
t0_p = tFWHM_p/2/sqrt(log(2));
C2_p = 1.3e-23;
l1 = l0;
l2 = 810e-9;
f_p = c/l1;
tshift = 0;
w0 = 2*pi*f_p;
P_p = 0.94*(J_p/tFWHM_p);
Epump = gaussianPulse(P_p,C2_p,t0_p,f_p, tshift, t, f, f0);
singlePlot(Epump, t, lbd, lambda_low, lambda_high, 'log')
spec = fftshift(fft(Epump)*dt)/1e-12;
% phi_f = angle(spec);
% phi_t = angle(Epump);
% 
fR = 0.18;
% 
gamma_SUP5 = 0.012;
% gamma_ESM10 = 0.0012;
alpha = 0;
% betas_SUP5 = beta_from_D('D_fibre_cesta_modif_1060.txt',l0);
betas_SUP5 = [1.71824211662369e-27, 6.46632805049671e-41,...
    -9.20155975986217e-56, 1.35083264938213e-70, 5.32762338120403e-85...
    , 5.02660485137516e-98, -7.51409193362055e-112, ...
    -2.04603524066773e-125, 3.67882131528204e-139,...
    6.22380068318912e-153, -1.27663429922201e-166];
% 
% betas_ESM10 = [3.27976403564214e-26, 3.09940609846371e-41, ...
%     -1.81813722695593e-57, 6.80603292588619e-70, 1.72134411055873e-83, ...
%     3.42174841826936e-97, 4.92290391762715e-111, 5.09548566219803e-125, ...
%     3.61458797820265e-139, 1.58533597626118e-153];
% 
% 
L_SUP5 = 10e-2;
h = L_SUP5/200;
tol = 1e-5;
Eout = rectPulse(0,Tspan, f0, t, f, f0);
E = Epump;
% Esave = zeros(100,res);
% for jk = 1:100
%     fprintf(1, ['Roundtrip no : ' num2str(jk) '\n']);
%     [Ecoupl, ~, ~] = propagationFibre(Eout, 1, 1, l0, l0, tol, t, f, lbd, 9,...
%                 [0], 0, 0, 'Coupler');
%     Efiltre = spectralFilter(Ecoupl, 950e-9, 40e-9, t, lbd, 'LPF');
%     [Efeedback, ~, ~] = propagationFibre(Efiltre, 210,210/20, l0, l2, 1e-4, t, f, lbd,...
%                 alpha, betas_ESM10, gamma_ESM10, 0, 'ESM10');
% %     singlePlot(Efeedback, t, lbd, lambda_low, lambda_high, 'linear');
%     E = Efeedback+Epump;
    [Eout, ~, ~] = propagationFibre(E, L_SUP5, h, l0, l1, tol, t, f, lbd,...
                alpha, betas_SUP5, gamma_SUP5, fR, 'SUP5');
    singlePlot(Eout, t, lbd, lambda_low, lambda_high, 'log');
%     Esave(jk,:) = Eout;
%     save('datas.mat','t','f','lbd','Esave')
% end