function [sig_a, sig_e] = VLMAcrossSections(lbd)

Ai_a = [215.7, 503, 238.8, 1951, 97.26, 20.57]*1e-27;
Bi_a = [902.7, 920, 971.6, 976.8, 1003, 1018]*1e-9;
Ci_a = [41.95, 30.15, 18.46, 4.44, 33.19, 14.59] *1e-9;
Ai_e = [1.03, 46.89, 45.19, 1961, 204.9, 294.11, 480.9, 70.81]*1e-27;
Bi_e = [851.7, 936.8, 970, 977.3, 980, 1018, 1027, 1070]*1e-9;
Ci_e = [1.23, 26.91, 23.25, 4.41, 15.47, 41.09, 21.93, 48.29] *1e-9;

sig_a = zeros(1,length(lbd));
sig_e = zeros(1,length(lbd));

for jk = 1:length(Ai_a)
    sig_a =sig_a+(Ai_a(jk).*exp(-((lbd-Bi_a(jk))./Ci_a(jk)).^2));
end

for jk = 1:length(Ai_e)
    sig_e = sig_e + (Ai_e(jk).*exp(-((lbd-Bi_e(jk))./Ci_e(jk)).^2));
end