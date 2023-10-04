% BME599 F23 | HW1 P2b(i)
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023

function Signals = p2b1()

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
phi_spoiler = ([1:N]/N-0.5 ) * 8*pi;

% Sweep off resonance frequencies
frequencies = -100:1:100;
nfreq = length(frequencies);

% To store the signal from each spin in the voxel
Signals = zeros(nfreq,1);

for ff=1:length(frequencies)
    dfreq = frequencies(ff);
%     disp(dfreq);

    M = zeros(3,N);

    for k = 1:N
	    [M1sig,M1] = gssignal(alpha,T1,T2,TE,TR,dfreq,phi_spoiler(k));
	    M(:,k)=M1;
    end
    
    Mss = mean(M')';
    Msig = Mss(1)+i*Mss(2);
    Signals(ff) = Msig;

end

magnSignals = abs(Signals);
phaseSignals = angle(Signals);

if 0 == 1
    f = figure('position',[495 473 954 331],'color','w'); 
    subplot(121);
    plot(frequencies, magnSignals);
    xlabel('Off-resonance frequencies (Hz)');
    ylabel('SS signal (magnitude)');
    title(['Gradient spoiled steady state simulation - N = ' num2str(N) ' spins']);
    ylim([0 1]);
    
    subplot(122);
    plot(frequencies, phaseSignals);
    xlabel('Off-resonance frequencies (Hz)');
    ylabel('SS signal (phase)');
    title(['Gradient spoiled steady state simulation - N = ' num2str(N) ' spins']);

end

