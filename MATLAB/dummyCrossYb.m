function [sig_a, sig_e] = dummyCrossYb(lbd)

Ai_a = [299.92*1e-27, 2314.71*1e-27, 75.54*1e-27, 211.61*1e-27, 496.62*1e-27, 5.12*1e-27];
Bi_a = [973.12*1e-9, 977.05*1e-9, 1014.99*1e-9, 908.36*1e-9, 917.03*1e-9, 1059.54*1e-9];
Ci_a = [32.84*1e-9, 7.24*1e-9, 37.69*1e-9, 92.3*1e-9, 48.67*1e-9, 40.2*1e-9];
Ai_e = [268.98*1e-27, 633.48*1e-27, 104.74*1e-27, 323.84*1e-27, 2299.03*1e-27];
Bi_e = [981.68*1e-9, 1025.96*1e-9, 975.32*1e-9, 1072.08*1e-9, 977.53*1e-9];
Ci_e = [27.16*1e-9, 42.01*1e-9, 86.16*1e-9, 44.96*1e-9, 7.25*1e-9];

sig_a = zeros(1,length(lbd));
sig_e = zeros(1,length(lbd));

for jk = 1:length(Ai_a)
    sig_a =sig_a+(Ai_a(jk).*exp(-((lbd-Bi_a(jk))./Ci_a(jk)).^2));
end

for jk = 1:length(Ai_e)
    sig_e = sig_e + (Ai_e(jk).*exp(-((lbd-Bi_e(jk))./Ci_e(jk)).^2));
end