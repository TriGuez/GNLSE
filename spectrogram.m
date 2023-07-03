function Sp = spectrogram(P,G,dt,N)

npt = length(P);
step_idx = round(npt/N);
N = ceil(npt/step_idx);
Sp = zeros(N,N);

P = P(:);
G = [zeros(npt/2-1,1);G(:);zeros(npt/2+1,1)];
idx_t = 1:1:npt;

kl = 1;

for idx_taux = 1:step_idx:npt
    Gprime = G(idx_t-(idx_taux-npt));
    E = abs(fftshift(fft(fftshift(P.*Gprime)))).^2.*dt^2;
    Sp(:,kl) = E(1:step_idx:end);
    kl = kl +1;
end