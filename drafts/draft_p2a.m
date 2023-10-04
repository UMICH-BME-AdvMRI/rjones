


% Bloch Equation Simulation, Excercise A-5b
% -----------------------------------------
% 

% dT = 1;		% 1ms delta-time.
% T = 1000;	% total duration
% N = ceil(T/dT)+1; % number of time steps.
% df = 10;	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.


TE = 2.5;
TR = TE*2;
alpha = deg2rad(60);

freqrange = -200:1:200;
nfreq = length(freqrange);

clear freqresp
freqresp.sig.complex = zeros(nfreq,1);
freqresp.ss = zeros(nfreq,3);

for ii = 1:length(freqrange)
    dfreq = freqrange(ii);

    [Msig,Mss] = sssignal(alpha,T1,T2,TE,TR,dfreq);
    freqresp.ss(ii,:) = Mss;
    freqresp.sig.complex(ii,1) = Msig;

end

freqresp.sig.magn = abs(freqresp.sig.complex);

figure; plot(freqresp.sig.magn);




%%
T1 = 300;	% ms.
T2 = 100;	% ms.
alpha_deg = 60;
alpha = deg2rad(alpha_deg);
frequencies = -200:1:200;

res1 = sim_ss_freq_response( alpha, T1, T2, 2.5, 5, frequencies );
res2 = sim_ss_freq_response( alpha, T1, T2, 5, 10, frequencies );
res3 = sim_ss_freq_response( alpha, T1, T2, 10, 20, frequencies );

% res1 = sim_ss_freq_response_x( alpha, T1, T2, 2.5, 5, frequencies );
% res2 = sim_ss_freq_response_x( alpha, T1, T2, 5, 10, frequencies );
% res3 = sim_ss_freq_response_x( alpha, T1, T2, 10, 20, frequencies );

f = figure('color','w'); % ,'position',[]);
hold on;
plot(frequencies,res1.sig.magn,'LineWidth',2);
plot(frequencies,res2.sig.magn,'LineWidth',2);
plot(frequencies,res3.sig.magn,'LineWidth',2);
xlim([min(frequencies) max(frequencies)])
xlabel('Frequency (Hz)');
ylabel('Steady state signal magnitude')
title(['bSSFP simulation' ...
    ' (\alpha=' sprintf('%d',alpha_deg) '^{\circ}, ' ...
    'T1=' sprintf('%d',T1) ', ' ...
    'T2=' sprintf('%d',T2) ')']);




