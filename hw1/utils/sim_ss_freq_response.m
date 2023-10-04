function res = sim_ss_freq_response( alpha, T1, T2, TE, TR, frequencies )
% function res = sim_ss_freq_response( alpha, T1, T2, TE, TR, frequencies )
% 
% Simulates steady state signal frequency response for given alpha (flip
% angle), T1, T2, TE and TR across specified off-resonances (frequencies)
%
% INPUTS:
%   alpha       = flip angle (radians)
%   T1          = T1 relaxation time (ms)
%   T2          = T2 relaxation time (ms)
%   TE          = echo time (ms)
%   TR          = repetition time (ms)
%   frequencies = off-resonance frequencies to simulate signal for (Hz)
%
% OUTPUTS:
%   res         = struct containing fields with steady state signal,
%                 magnetization components and simulation parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE NOTE:
%  Uses functions (e.g., sssignal()) from Brian Hargreaves' Bloch Simulation
%   website [ http://mrsrl.stanford.edu/~brian/bloch/ ]
%
% RJ 10-03-2023

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
res.sig.phase = angle(res.sig.complex);

res.params.alpha = alpha;
res.params.T1 = T1;
res.params.T2 = T2;
res.params.TE = TE;
res.params.TR = TR;
res.params.frequencies = frequencies;