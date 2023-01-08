function TE = linearProp(E, f, l0, lc, betas, L)

c = 299792458;
f0 = c./l0;
fc = c./lc;
f_cent = (f+f0)-fc;

wshift = fftshift(2*pi*f_cent);

B = betas(1)./2*wshift.^2;
for jl = 2:length(betas)
    B = B + betas(jl)/factorial(jl+1) * wshift.^(jl+1);
end

TE = ifft(fft(E).*exp(-1i*B*L));