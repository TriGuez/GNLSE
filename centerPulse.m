function Ecenter = centerPulse(E, t)

idx = find(abs(E).^2 == max(abs(E).^2));
N = length(t)/2;
idx = (N-idx);

Ecenter = circshift(E,idx);