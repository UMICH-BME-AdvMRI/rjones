% BME599 F23 | HW1 P2b(i)
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023
%
% Script for HW1 Problem 2 part b(i)
%

function Signals = p2b1_gssignal(makePlots)
% Signals = p2b1(makePlots)
%
% runs simulations for hw1 p2 part(i)
% Set makeplots=1 to create + save plots of results
% [made this a function so it could be called(+return Signals) elsewhere]
%

if nargin<1, makePlots=0; end

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

% Spoiler "phase twists" to apply to spins in a voxel
%  ( total phase twist = 8pi )
phi_spoiler = ([1:N]/N-0.5 ) * 8*pi;

% To sweep through grad spoiling frequencies
frequencies = -100:1:100;
nfreq = length(frequencies);

% To store the s.s. signals (for each off-resonance freq)
Signals = zeros(nfreq,1);
MSignals = zeros(nfreq,3);

for ff=1:length(frequencies)
    dfreq = frequencies(ff);
%     fprintf(' orr-resonance-freq=%d\n',dfreq);

    % mx,my,mz for all N spins in voxel
    M = zeros(3,N);
    % simulation for each spin in voxel
    for k = 1:N
	    [M1sig,M1] = gssignal(alpha,T1,T2,TE,TR,dfreq,phi_spoiler(k));
	    M(:,k)=M1;
    end
    % average across spins in voxel
    Mss = mean(M')';
    % calculate M_xy transverse magnetization
    Msig = Mss(1) + 1i*Mss(2);
    % store result
    Signals(ff) = Msig;
    MSignals(ff,:) = Mss;
end

magnSignals = abs(Signals);
phaseSignals = angle(Signals);

if makePlots
  % Plot the signal magnitude and phase as a function of off-resonance freq
    f = figure('position',[495 88 505 716],'color','w'); 
    subplot(211);
    plot(frequencies, magnSignals,'LineWidth',2);
    xlabel('Off-resonance frequency (Hz)');
    ylabel('SS signal (magnitude)');
    title(['Gradient spoiled steady state simulation - N = ' num2str(N) ' spins']);
    ylim([0 1]);
    set(gca,'FontSize',15);
    annotation(f,'textarrow',[0.552475247524753 0.520792079207921],...
    [0.677770949720671 0.613128491620112],'String',...
    {['ss signal = ' sprintf('%g',magnSignals(101))],'(for all off-resonance frequencies)'},'FontSize',15);

    
    subplot(212);
    plot(frequencies, phaseSignals,'LineWidth',2);
    xlabel('Off-resonance frequency (Hz)');
    ylabel('SS signal (phase)');
    title(['Gradient spoiled steady state simulation - N = ' num2str(N) ' spins']);
    set(gca,'FontSize',15);

%     fout = 'plots/p2b-i-SignalMagn+Phase.png';
%     print(f,fout,'-dpng');


  % Plot the magnetization components M_{x,y,z} w.r.t. off-res-freq
    f = figure('position',[270 407 539 404],'color','w'); 
    hold on;
    plot(frequencies, MSignals(:,1),'LineWidth',2);
    plot(frequencies, MSignals(:,2),'LineWidth',2);
    plot(frequencies, MSignals(:,3),'LineWidth',2);
    xlabel('Off-resonance frequency (Hz)');
    ylabel('Residual magnetization');
    lgn = legend({'M_x','M_y','M_z'});
    title(['Gradient spoiled steady state simulation - N = ' num2str(N) ' spins']);
%     ylim([0 1]);
    set(gca,'FontSize',15);
    
%     fout = 'plots/p2b-i-Mxyz.png';
%     print(f,fout,'-dpng');


end

