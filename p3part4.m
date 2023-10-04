% -----------------------------------------
% BME599 F23 | HW1 P3(4) - Slice Profile Simulation
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023
% -----------------------------------------

clear
close all
% clc


% Gamma bar [Hz/T] (gyromagnetic ratio)
gyro = 42.58;  % MHz/T

% RF pulse duration
tau_rf = 2 * 1e-3; % [s]
% (discrete) time points
dt = 5e-3 * 1e-3; % time step size [s]
% (discrete) time vector
t = 0:dt:tau_rf; % [s]

%%%% Create RF pulse waveform:

% Set time bandwidth product
TBW = 8;
BW = TBW/tau_rf; % [Hz]

% flip angle
alpha = 90; % degrees

% -rf shape- sinc with # of zero-croxxings == TBW
% -rf amplitude- make such that integral of rf waveform gives desired flip
%                   angle (here we want 90 deg filp)
t_rf = 1e3*(t - (tau_rf/2));
rf = sinc(TBW/2 * t_rf);  % sinc(t) = sin(pi*t)/(pi*t)
rf_amp = 1e3 * alpha/(360*gyro*sum(rf));
rf90 = rf_amp * rf;
rf_orig90 = rf;

% Spatial profile
dz = 0.1; % mm
z = -20:dz:20;  % mm

% Sequence params
T1 = 1000;      % ms
T2 = 100;       % ms

% G_sliceselect
grad_amp = 18.788; % mT/m
grad_amp = grad_amp * 1e-3 * 1e4 * 1e-2; % Gauss/cm
grad = grad_amp * ones(size(t));

%gradient rise time
slice_thickness = 5*1e-3; % [m]
slew_rate = 180; % mT/m/ms
tau_grad_rise = ((BW/(gyro*1e6 * slice_thickness)))/slew_rate; % [s] 
slew_rate = slew_rate * 1e-3 * 1e4 * 1e-2 * 1e3; % Gauss/cm/s

nrise = round(tau_grad_rise/dt);
grad_rise_up = linspace(0,1,nrise)*grad_amp;
grad_rise_down = flip(grad_rise_up,2);

t_buffer_s = -nrise*dt:dt:0-dt;
t_buffer_e = tau_rf+dt:dt:tau_rf+(nrise*dt);

% add in gradient ramping
grad_full = cat(2,grad_rise_up,grad,grad_rise_down);
t_full = cat(2,t_buffer_s,t,t_buffer_e);
t_full = t_full-min(t_full);
rf90_full = cat(2,0*grad_rise_up,rf90,0*grad_rise_down);

% %%% Add rephasing gradient (no grad ramping)
% grad_re = [grad -0.5*grad];
% rf90_re = [rf90 zeros(size(t))];
% t_re    = [t t+t(end)+dt];

%%% Add rephasing gradient (with grad ramping)
stopIndex = length(grad_full);
grad_re = [grad_full -0.5*grad_full];
t_re    = [t_full t_full+t_full(end)+dt];
rf90_re = [rf90_full zeros(size(t_full))];

% Zero off-resonance
df = 0;

% Simulation at time between slice-select grad and slice-rephasing grad
% ( i.e., simulation without slice-rephasing gradient included)
[msig, m] = sliceprofile_part3(rf90_re,grad_re,t_re,T1,T2,z,df,stopIndex);

% Simulation at time AFTER slice-rephasing gradient
[msig_re, m_re] = sliceprofile(rf90_re,grad_re,t_re,T1,T2,z,df);


%% Plot signal magnitude and phase

f = figure('color','w','position', [476 240 595 626]);
subplot(211);
hold on;
plot(z,abs(msig),'LineWidth',2);
plot(z,abs(msig_re),'LineWidth',2);
title('M_{xy} Magnitude Slice Profile (flip=90deg, df=0Hz)');
legend({'Without slice rephasing','With slice rephasing'});
xlabel('Position (mm)');
ylabel('Signal Magnitude');
ylim([0 1]);
grid on;
set(gca,'FontSize',15);
subplot(212);
hold on;
plot(z,angle(msig),'LineWidth',2);
plot(z,angle(msig_re),'LineWidth',2);
title('M_{xy} Phase Slice Profile (flip=90deg, df=0Hz)');
legend({'Without slice rephasing','With slice rephasing'});
xlabel('Position (mm)');
ylabel('Phase (rad)');
ylim([-pi pi]);
grid on;
set(gca,'FontSize',15);


fout = 'plots/p3-part4-slice-rephasing-SignalMagn+Phase-comparison.png';
print(f,fout,'-dpng');



%% Plot Mx, my, mz 

f = figure('color','w','position',[438 175 525 677]);
subplot(211);
hold on;
plot(z, m(1,:), 'LineWidth', 2);
plot(z, m(2,:), 'LineWidth', 2);
plot(z, m(3,:), 'LineWidth', 2);
plot(z, zeros(size(z)),'-k');
ylim([-1 1]);
set(gca,'FontSize',14);
xlabel('Position (mm)');
ylabel('Residual Magnetization');
lgn = legend({'M_x','M_y','M_z'}); %,'Position',[0.707662271805274 0.727282658293597 0.114604462474645 0.0774032459425719]);
title({'P3.4 | M profiles | No slice rephsaing'}); % ...
%     ['TBW=' num2str(TBW) ', \tau_{RF}=' num2str(tau_rf)]});
subplot(212);
hold on;
plot(z, m_re(1,:), 'LineWidth', 2);
plot(z, m_re(2,:), 'LineWidth', 2);
plot(z, m_re(3,:), 'LineWidth', 2);
plot(z, zeros(size(z)),'-k');
ylim([-1 1]);
set(gca,'FontSize',14);
xlabel('Position (mm)');
ylabel('Residual Magnetization');
lgn = legend({'M_x','M_y','M_z'}); %,'Position',[0.707662271805274 0.727282658293597 0.114604462474645 0.0774032459425719]);
title({'P3.4 | M profiles | With slice rephsaing'}); % ...

fout = 'plots/p3-part4-slice-rephasing-Mxyz-comparison.png';
print(f,fout,'-dpng');






