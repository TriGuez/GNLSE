function E = rectPulse(P, tFWHM, f1, t, f, f0)

%  This function computes the complex enveloppe of a rectangular optical pulse
%  centered at frequency f1 and optical noise of one photon per spectral node
%  INPUTS : 
%        P : Peak power of the pulse [W]
%        tFWHM : Half temporal width of the pulse [s]
%        f1 : Central frequency of the pulse [Hz]
%        t : Time vector of the simulation [s]
%        f : Frequency vector of the pulse [Hz]
%        f0 : Central frequency of the simulation [Hz]
%  OUTPUTS :
%        E : Complex enveloppe of the optical pulse 

dt = t(2)-t(1);
df = f(2)-f(1);

res = length(t);
h_bar = 6.6207004e-34;
noise = sqrt(h_bar*fftshift(f+f0)/df*(res-1)/res).*exp(-1i*rand(1,res)*2*pi);
noise(isnan(noise)) = 0+0*1i;
spectralPhase = zeros(1,res);
E = sqrt(P).*((t >-tFWHM) & (t<tFWHM)).*exp(-1i*2*pi*(f0-f1)*t);
E = ifft(fft(E)*dt+noise)./dt;