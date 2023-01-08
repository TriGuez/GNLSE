function E = gaussianPulse(P,C2,t0,f1,t_shift,t,f, f0)

% This function computes a gaussian optical pulse 
% INPUTS : 
%       P : Peak power [W]
%       C2 : 2nd order spectral phase [s**2]
%       t0 : 1/e half pulse duration [s]
%       f1 : Optical center frequency [Hz]
%       t_shift : Time domain shift [s]
%       t : Simulation time vector [s]
%       f : Simulation frequency vector [Hz]
%       f0 : Simulation central frequency [Hz]
% OUTPUTS : 
%       E : Corresponding gaussian optical pulse


dt = t(2)-t(1);
df = f(2)-f(1);
res = length(f);
c = 299792458;
h = 6.62607004e-34;
wshift = fftshift(2*pi*f);

noise = sqrt(h*fftshift(f+f0)/df*(res-1)/res).*exp(-1i*rand(1,res)*2*pi);
% noise = zeros(1,length(noise));
spectral_phase = C2/2*wshift.^2;
E = sqrt(P).*exp(-0.5*(((t-t_shift)/t0).^2)).*exp(-1i*2*pi*(f0-f1).*t);
E = ifft(fft(E).*dt.*exp(-1i*spectral_phase)+noise)/dt;