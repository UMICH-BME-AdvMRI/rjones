function res = sim_ss_freq_response( alpha, T1, T2, TE, TR, frequencies )

nfreq = length(frequencies);
res.sig.complex = zeros(nfreq,1);
res.ss = zeros(nfreq,3);
for ii = 1:nfreq
    dfreq = frequencies(ii);
    [Msig,Mss] = sssignal(alpha,T1,T2,TE,TR,dfreq);
    res.ss(ii,:) = Mss;
    res.sig.complex(ii,1) = Msig;
end
res.sig.magn = abs(res.sig.complex);

res.params.alpha = alpha;
res.params.T1 = T1;
res.params.T2 = T2;
res.params.TE = TE;
res.params.TR = TR;
res.params.frequencies = frequencies;