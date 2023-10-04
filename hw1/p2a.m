% BME599 F23 | HW1 P2A(1-3)
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023


%% Run this to add the necessary directories+files to search path

% (( may need to adjust path
% addpath(genpath('/Users/robertjones/Desktop/F23/599/hw/hw1/blochsim_matlab'));
addpath(genpath(pwd));


%% Set parameters

% Pick some values for T1, T2 (not specified in problem)
T1 = 1000;	% [ms]
T2 = 100;	% [ms]

% Set specified flip angle
alpha = 60; % [deg]

% Pick some (off-resonance) frequency range to evaluate ss signal response
frequencies = -200:1:200; % [Hz]

% Set TRs & TEs [given in problem statement]
TRs = [5, 10, 20];  % milliseconds [ms]
TEs = TRs/2;        % [ms]


%% Run bloch simulation

% Simulate the s.s. freq. response for TR/TE values given in P2A(1-3)
res1 = sim_ss_freq_response( deg2rad(alpha), T1, T2, TEs(1), TRs(1), frequencies );
res2 = sim_ss_freq_response( deg2rad(alpha), T1, T2, TEs(2), TRs(2), frequencies );
res3 = sim_ss_freq_response( deg2rad(alpha), T1, T2, TEs(3), TRs(3), frequencies );


%% Plots results

% Plot frequency responses for parts (1-3)
f = figure('color','w'); % ,'position',[]);
hold on;
plot(frequencies,res1.sig.magn,'LineWidth',2); %te/tr=2.5/5 ms
plot(frequencies,res2.sig.magn,'LineWidth',2); %te/tr=5/10 ms
plot(frequencies,res3.sig.magn,'LineWidth',2); %te/tr=10/20 ms
xlim([min(frequencies) max(frequencies)]);
xlabel('Off-resonance frequency (Hz)');
ylabel('Steady state signal magnitude')
title(['bSSFP simulation' ...
    ' (\alpha=' sprintf('%d',alpha) '^{\circ}, ' ...
    'T1=' sprintf('%d',T1) 'ms, ' ...
    'T2=' sprintf('%d',T2) 'ms)']);
lgn = legend({'TR/TE = 5/2.5 ms','TR/TE = 10/5 ms','TR/TE = 20/10 ms'},'location','best');
set(gca,'FontSize',14);
drawnow; shg;

fout = ['plots/p2a-ss-Magn-freqresponse_T1,' num2str(T1) '-T2,' num2str(T2) '.png'];
print(f,fout,'-dpng');
close(f);



