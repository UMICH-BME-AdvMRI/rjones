% BME599 F23 | HW1 P2b(iii)
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023
%
% Script for HW1 Problem 2 part b(iii)
%

clear
close all
% clc

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
ncycles = 4; % 2*pi rad [e.g., ncycles=4 == 8pi dephasing]

% Set # of RF excitations (kinda like "train length"?)
Nex = 100;
% Set # of frequencies
Nf = 100;

% Off-resonance frequency
dfreq = 0;

% RF pulse phases to test [every 1degree from 0 to 180 deg]
thetas = deg2rad(0:1:180);

% To store the signal at each RF pulse phase
S = zeros(length(thetas),1);
% To store the signal from each excitation  at each RF pulse phase
S_AllEx = zeros(length(thetas),Nex);
% TO store magnetization components  at each RF pulse phase
MS = zeros(length(thetas),3,Nex);

% Loop through RF pulse phases
legendstring = {};
for ind=1:length(thetas)
    rfphase = thetas(ind);
%     [Msig,Mss] = gresignal(alpha,T1,T2,TE,TR,dfreq,ncycles,N);
%     [Msig,Mss] = gssignal_rfspoil(alpha,T1,T2,TE,TR,dfreq,ncycles,rfphase);
%     [Msig,Mss] = gresignal_rfspoil(alpha,T1,T2,TE,TR,dfreq,ncycles,N,rfphase);
    [Msig,Mss,Msigs]=spgrsignal(alpha,T1,T2,TE,TR,dfreq,Nex,rfphase,Nf,ncycles);

    S(ind) = Msig;
    S_AllEx(ind,:) = Msigs;
    MS(ind,:,:) = Mss;

    legendstring{ind} = num2str(rad2deg(rfphase));
end

magnSignals = abs(S);
phaseSignals = angle(S);

% Get the "perfect"/target signal from gradient spoiled simulation (i.e.,
%       the average across the entire voxel)
Signals = p2b1();
grspsignal = mean(abs(Signals));

f=figure('color','w'); %,'position',[]);
hold on;
plot(rad2deg(thetas),abs(S),'-','LineWidth',2,'Marker','o'); 
plot(rad2deg(thetas),repmat(grspsignal,length(thetas),1),'--m','LineWidth',2);
xlabel('RF pulse phase [deg]');
ylabel('SS signal magnitude');
title('P2b(iii) | Grad+RF spoiled SS simulations');

shg


f = figure('position',[495 473 954 331],'color','w'); 
hold on;
plot(rad2deg(thetas), MS(:,1,end) ,'-o','LineWidth',2,'Marker','o');
plot(rad2deg(thetas), MS(:,2,end) ,'-o','LineWidth',2,'Marker','o');
plot(rad2deg(thetas), MS(:,3,end) ,'-o','LineWidth',2,'Marker','o');
plot(rad2deg(thetas), abs(MS(:,1,end)+1i*MS(:,2,end)) ,'-o','LineWidth',2,'Marker','o');

xlabel('RF pulse phase (rad)');
ylabel('SS residual magnitude');
title('P2b(iii) | Grad+RF spoiled SS simulations');


% 
% f = figure('position',[495 473 954 331],'color','w'); 
% hold on;
% plot(rad2deg(thetas(1:10:end)), abs(S_AllEx(1:10:end,:)).' ,'-','LineWidth',2);
% 
% plot(rad2deg(thetas), MS(:,2,end) ,'-o','LineWidth',2,'Marker','o');
% plot(rad2deg(thetas), MS(:,3,end) ,'-o','LineWidth',2,'Marker','o');
% plot(rad2deg(thetas), abs(MS(:,1,end)+1i*MS(:,2,end)) ,'-o','LineWidth',2,'Marker','o');
% 
% xlabel('RF pulse phase (rad)');
% ylabel('SS residual magnitude');
% title('P2b(iii) | Grad+RF spoiled SS simulations');
% 
% f = figure('position',[495 473 954 331],'color','w'); 
% plot(rad2deg(thetas(117:119)), abs(S_AllEx(117:119,:)) ,'-','LineWidth',2);


