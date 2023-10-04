% BME599 F23 | HW1 P2A(1-3)
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023

%%%%%%%%    Problem 2(a):       %%%%%%%%
% 2a. Simulate the steady-state frequency response of a bSSFP sequence with
% a flip angle of 60° using Bloch equation simulations. Show results for
% (1) TR=5ms and TE=2.5ms, (2) TR=10ms and TE=5ms, and (3) TR=20ms and TE=10ms.

%% Set parameters

% Pick some values for T1, T2 (not specified in problem)
T1 = 300;	% ms.
T2 = 100;	% ms.

% Set specified flip angle
alpha = deg2rad(60); % deg
% Get rotation matrix for RF excitation flip (rotate about x or y)
%  [this tips magnetization into the transverse plane by angle alpha]
R_tip = xrotation(alpha);

% Pick some frequency range to eval ss response
frequencies = -200:1:200; % Hz

% TR values to use
TRs = [5, 10, 20];  % ms

% TE values to use
TEs = TRs/2;

% to store the ss signals
ntrials = length(TRs);
nfreq = length(frequencies);
results = zeros(nfreq, ntrials);

for jj=1:ntrials
    TR = TRs(jj);
    TE = TRs(jj);
    for kk=1:nfreq
        df = frequencies(kk);  % off-resonance frequency
        
        phi = 2*pi*df*t/1000;	% Resonant precession, radians.
        E1 = exp(-t/t1);	
        E2 = exp(-t/t2);
        
        Afp = [E2 0 0;0 E2 0;0 0 E1]*zrot(phi);
        Bfp = [0 0 1-E1]';
    end
end



% nfreq = length(frequencies);
% res.sig.complex = zeros(nfreq,1);
% res.ss = zeros(nfreq,3);
% for ii = 1:nfreq
%     dfreq = frequencies(ii);
%     [Msig,Mss] = sssignal(alpha,T1,T2,TE,TR,dfreq);
%     res.ss(ii,:) = Mss;
%     res.sig.complex(ii,1) = Msig;
% end
% res.sig.magn = abs(res.sig.complex);


[Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);
[Ate,Bte] = freeprecess(TE,T1,T2,dfreq);
Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
Msig = Mss(1)+i*Mss(2);
% Let 	M1 be the magnetization just before the tip.
%	M2 be just after the tip.
%	M3 be at TE.
%
% then
%	M2 = Rflip * M1
%	M3 = Ate * M2 + Bte
%	M1 = Atr * M3 + Btr
%
% Solve for M3...
%
%	M3 = Ate*Rflip*Atr*M3 + (Ate*Rflip*Btr+Bte)


phi = 2*pi*df*T/1000;	% Resonant precession, radians.
E1 = exp(-T/T1);	
E2 = exp(-T/T2);
Afp = [E2 0 0;0 E2 0;0 0 E1]*zrot(phi);
Bfp = [0 0 1-E1]';




%% OLD:

% simulate ss response for TR/TE values for parts (1-3)
res1 = sim_ss_freq_response( deg2rad(alpha), T1, T2, 2.5, 5, frequencies );
res2 = sim_ss_freq_response( deg2rad(alpha), T1, T2, 5, 10, frequencies );
res3 = sim_ss_freq_response( deg2rad(alpha), T1, T2, 10, 20, frequencies );

% res1 = sim_ss_freq_response_x( alpha, T1, T2, 2.5, 5, frequencies );
% res2 = sim_ss_freq_response_x( alpha, T1, T2, 5, 10, frequencies );
% res3 = sim_ss_freq_response_x( alpha, T1, T2, 10, 20, frequencies );

% Plot frequency responses for parts (1-3)
f = figure('color','w'); % ,'position',[]);
hold on;
plot(frequencies,res1.sig.magn,'LineWidth',2); %te/tr=2.5/5 ms
plot(frequencies,res2.sig.magn,'LineWidth',2); %te/tr=5/10 ms
plot(frequencies,res3.sig.magn,'LineWidth',2); %te/tr=10/20 ms
xlim([min(frequencies) max(frequencies)]);
xlabel('Frequency (Hz)');
ylabel('Steady state signal magnitude')
title(['bSSFP simulation' ...
    ' (\alpha=' sprintf('%d',alpha) '^{\circ}, ' ...
    'T1=' sprintf('%d',T1) ', ' ...
    'T2=' sprintf('%d',T2) ')']);
lgn = legend({'TR/TE = 5/2.5 ms','TR/TE = 10/5 ms','TR/TE = 20/10 ms'},'location','best');
set(gca,'FontSize',14);
drawnow; shg;
print(f,'figures/p2a.png','-dpng');
close(f);



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



