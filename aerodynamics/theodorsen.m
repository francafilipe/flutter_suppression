function C = theodorsen(reduced_freq)

C = zeros(length(reduced_freq),1);
for k=1:length(reduced_freq)

    if reduced_freq(k)==0
        reduced_freq(k) = 1e-6;
    end

    K1 = besselk(1,1i*reduced_freq(k));
    K0 = besselk(0,1i*reduced_freq(k));

    C(k) = K1/(K0+K1);
end

end