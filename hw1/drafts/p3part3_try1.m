% -----------------------------------------
% BME599 F23 | HW1 P3(3) - Slice Profile Simulation
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
tau_rf = 2*1e-3; % [s]
% (discrete) time points
dt = 5e-6; % time step size [s]
% (discrete) time vector
t = 0:dt:tau_rf; % [s]


%%%% Create RF pulse waveform:

% Set time bandwidth product
TBW = 8;


%%% Generate 10deg flip angle RF pulse - 
% flip angle
alpha = 90; % degrees

% -rf shape- sinc with # of zero-croxxings == TBW
% -rf amplitude- make such that integral of rf waveform gives desired flip
%                   angle (here we want 90 deg filp)
t_rf = 1e3*(t - (tau_rf/2));
rf = sinc(TBW/2 * t_rf);  % sinc(t) = sin(pi*t)/(pi*t)
rf_amp = 1e3 * alpha/(360*gyro*sum(rf));
rf90 = rf_amp * rf;


%%% Generate 10deg flip angle RF pulse - 
% flip angle
alpha = 30; % degrees

% -rf shape- sinc with # of zero-croxxings == TBW
% -rf amplitude- make such that integral of rf waveform gives desired flip
%                   angle (here we want 90 deg filp)
t_rf = 1e3*(t - (tau_rf/2));
rf = sinc(TBW/2 * t_rf);  % sinc(t) = sin(pi*t)/(pi*t)
rf_amp = 1e3 * alpha/(360*gyro*sum(rf));
rf30 = rf_amp * rf;


%%% Generate 10deg flip angle RF pulse - 
% flip angle
alpha = 10; % degrees

% -rf shape- sinc with # of zero-croxxings == TBW
% -rf amplitude- make such that integral of rf waveform gives desired flip
%                   angle (here we want 90 deg filp)
t_rf = 1e3*(t - (tau_rf/2));
rf = sinc(TBW/2 * t_rf);  % sinc(t) = sin(pi*t)/(pi*t)
rf_amp = 1e3 * alpha/(360*gyro*sum(rf));
rf10 = rf_amp * rf;


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

%%% Add rephasing gradient
grad_re = [grad -0.5*grad];
rf90_re = [rf90 zeros(size(t))];
rf30_re = [rf30 zeros(size(t))];
rf10_re = [rf10 zeros(size(t))];
t_re    = [t t+t(end)+dt];

%%% off-resonance
df = 0; % Hz

%%% bloch simulation
[msig90,m90] = sliceprofile(rf90_re, grad_re, t_re, ...
    T1, T2, z, df);

%%% bloch simulation
[msig30,m30] = sliceprofile(rf30_re, grad_re, t_re, ...
    T1, T2, z, df);

%%% bloch simulation
[msig10,m10] = sliceprofile(rf10_re, grad_re, t_re, ...
    T1, T2, z, df);

% figure;
% subplot(211);
% plot(z,abs(msig30),'LineWidth',2);
% xlabel('Position (mm)');
% ylabel('Signal Magnitude');
% grid on;
% title('Magnitude Slice Profile - RF flip angle = 30deg');
% subplot(212);
% plot(z,abs(msig10),'LineWidth',2);
% xlabel('Position (mm)');
% ylabel('Signal Magnitude');
% grid on;
% title('Magnitude Slice Profile - RF flip angle = 10deg');


f = figure('color','w','position',[476 498 448 368]);
hold on;
plot(z,abs(msig90),'LineWidth',2);
plot(z,abs(msig30),'LineWidth',2);
plot(z,abs(msig10),'LineWidth',2);
title('Magnitude Slice Profile');
legend({'Flip angle = 90deg','Flip angle = 30deg','Flip angle = 10deg'});
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
set(gca,'FontSize',13);

fout = 'plots/p3-part3-flip-angle-comparison.png';
print(f,fout,'-dpng');


% 
% %%% make plots
% % make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, plotTag );
% make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, df )
% 
% %%% off-resonance
% df = 200; % Hz
% 
% %%% bloch simulation
% [msig,m] = sliceprofile(rf90_re, grad_re, t_re, ...
%     T1, T2, z, df);
% 
% %%% make plots
% % make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, plotTag );
% make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, df );
% 
% 
% 
% % % Plot RF waveform
% % f = figure('color','w'); %,'position',[]); 
% % hold on;
% % plot(t_rf, rf, 'LineWidth', 2);
% % plot(t_rf, zeros(size(t_rf)),'-k');
% % set(gca,'FontSize',13);
% % xlabel('Time (ms)');
% % ylabel('RF amplitude (\muT)');
% % 
% % title({'P3.1 | RF pulse waveform' ...
% %     ['TBW=' num2str(TBW) ', \tau_{RF}=' num2str(tau_rf)]});
% 
% 
% 
% figure;
% 
% subplot(3,2,1);
% plot(z,abs(msig));
% xlabel('Position (mm)');
% ylabel('Signal Magnitude');
% grid on;
% title('Magnitude Profile');
% 
% subplot(3,2,3);
% plot(z,angle(msig));
% xlabel('Position (mm)');
% ylabel('Signal Phase (rad)');
% grid on;
% title('Phase Profile');
% 
% subplot(3,2,5);
% plot(z,m(3,:));
% xlabel('Position (mm)');
% ylabel('Residual M_z');
% grid on;
% title('M_z Profile');
% 
% subplot(3,2,2);
% plot(t,rf90);
% xlabel('Time (s)');
% ylabel('RF (G)');
% grid on;
% title('RF vs Time');
% 
% subplot(3,2,4);
% plot(t,grad);
% xlabel('Time (s)');
% ylabel('Gradient (G/cm)');
% grid on;
% title('Gradient vs Time');
% 
% drawnow;
% 
% 
% figure;
% 
% subplot(1,3,1);
% plot(z,m(1,:));
% xlabel('Position (mm)');
% ylabel('M_x');
% grid on;
% title('M_x Profile');
% 
% subplot(1,3,2);
% plot(z,m(2,:));
% xlabel('Position (mm)');
% ylabel('M_y');
% grid on;
% title('M_y Profile');
% 
% subplot(1,3,3);
% plot(z,m(3,:));
% xlabel('Position (mm)');
% ylabel('M_z');
% grid on;
% title('M_z Profile');
% 
% drawnow;
% 
% 
% % % Plot RF waveform
% % f = figure('color','w'); %,'position',[]); 
% % hold on;
% % plot(t_rf, rf, 'LineWidth', 2);
% % plot(t_rf, zeros(size(t_rf)),'-k');
% % 
% % % Plot RF waveform
% % f = figure('color','w'); %,'position',[]); 
% % hold on;
% % plot(t, rf, 'LineWidth', 2);
% % plot(t, zeros(size(t_rf)),'-k');
% % yyaxis right; plot(t, grad_ss);

