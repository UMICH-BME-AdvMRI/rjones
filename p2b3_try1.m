% BME599 F23 | HW1 P2b(iii)
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023


%% Set parameters

% Number of spins in single voxel to simulate
N = 100;

% Flip angle
alpha = deg2rad(10); %rad

% T1/T2
T1 = 1000;  % ms
T2 = 100;   % ms

 % TR/TE
TR = 10;    % ms
TE = 5;     % ms

% Phase twists to apply to spins
%  ( total phase twist = 8pi )
ncycles = 4;

Nex = 100;
Nf = 100;

% Sweep off resonance frequencies
dfreq = 0;

% RF pulse phases to test
thetas = deg2rad(0:1:180);

% To store the signal from each spin in the voxel
S = zeros(length(thetas),1);
allS = zeros(length(thetas),Nex);


legendstring = {};
for ind=1:length(thetas)
    rfphase = thetas(ind);
%     [Msig,Mss] = gresignal(alpha,T1,T2,TE,TR,dfreq,ncycles,N);
%     [Msig,Mss] = gssignal_rfspoil(alpha,T1,T2,TE,TR,dfreq,ncycles,rfphase);
%     [Msig,Mss] = gresignal_rfspoil(alpha,T1,T2,TE,TR,dfreq,ncycles,N,rfphase);
    [Msig,Mss,Msigs]=spgrsignal(alpha,T1,T2,TE,TR,dfreq,Nex,rfphase,Nf,ncycles);


    S(ind) = Msig;
    allS(ind,:) = Msigs;

    legendstring{ind} = num2str(rad2deg(rfphase));
end

magnSignals = abs(S);
phaseSignals = angle(S);

Signals = p2b1();
grspsignal = mean(abs(Signals));

figure; 
hold on;
plot(rad2deg(thetas),abs(S)); 
plot(rad2deg(thetas),repmat(grspsignal,length(thetas),1),'--c')
shg


% f = figure('position',[495 473 954 331],'color','w'); 
% subplot(121);
% plot(thetas, magnSignals ,'-o');
% xlabel('RF pulse phase (rad)');
% ylabel('SS signal (magnitude)');
% title(['Gradient spoiled steady state simulation - N = ' num2str(N) ' spins']);
% % ylim([0 1]);
% 
% subplot(122);
% plot(frequencies, phaseSignals);
% xlabel('Off-resonance frequencies (Hz)');
% ylabel('SS signal (phase)');
% title(['Gradient spoiled steady state simulation - N = ' num2str(N) ' spins']);
% 


