% -----------------------------------------
% BME599 F23 | HW1 P3(1) - Slice Profile Simulation
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023
% -----------------------------------------

%%% Define known parameters
% Gamma bar [Hz/T] (gyromagnetic ratio)
gyro = 42.58*1e6;  % Hz/T
% RF pulse duration
tau_rf = 2; % ms
% Set time bandwidth product
TBW = 8; % a.u.
% flip angle
alpha = pi/2; % degrees
% T1/T2
T1=1000; % ms
T2=100; % ms
% G_ss amplitude
G_ss = 18.788; % mT/m

% (discrete) time points
dt = 0.001; % time step size [ms]
t = 0:dt:tau_rf; % ms
t_rf = t - tau_rf/2; % ms
t_seconds = t*1e-3;
nt = length(t);

BW = TBW/tau_rf;

%%%% Create RF pulse waveform:

% -rf shape- sinc with # of zero-croxxings == TBW
% -rf amplitude- make such that integral of rf waveform gives desired flip
%                   angle (here we want 90 deg filp)
rf = sinc(BW * t_rf);  % sinc(t) = sin(pi*t)/(pi*t)
rf_amp = 1e6 * (alpha / (2*pi * gyro * tau_rf*1e-3));  % microT
rf = rf_amp * (rf/sum(rf));
rf_ut = rf;
rf_gauss = rf_ut*1e1;

df = 0; % Hz
dx = 0.1;
xmax = 10;
x = (-xmax:dx:xmax); % mm
grad_amp = 18.788 * 1e-3 * 1e4 * 1e-2; % Gauss/cm
grad = grad_amp * ones(size(t));
grad_g_cm = grad;

%%% To add rephasing gradient in z/sliceselect dir: 
% grad = [grad -grad/2];
% t = [t t+t(end)];
% rf = [rf 0*rf];	


[msig,m] = sliceprofile(rf_gauss, grad_g_cm, t_seconds, T1, T2, x, df);





% % Plot RF waveform
% f = figure('color','w'); %,'position',[]); 
% hold on;
% plot(t_rf, rf, 'LineWidth', 2);
% plot(t_rf, zeros(size(t_rf)),'-k');
% set(gca,'FontSize',13);
% xlabel('Time (ms)');
% ylabel('RF amplitude (\muT)');
% 
% title({'P3.1 | RF pulse waveform' ...
%     ['TBW=' num2str(TBW) ', \tau_{RF}=' num2str(tau_rf)]});



figure;

subplot(3,2,1);
plot(x,abs(msig));
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
title('Magnitude Profile');

subplot(3,2,3);
plot(x,angle(msig));
xlabel('Position (mm)');
ylabel('Signal Phase (rad)');
grid on;
title('Phase Profile');

subplot(3,2,5);
plot(x,m(3,:));
xlabel('Position (mm)');
ylabel('Residual M_z');
grid on;
title('M_z Profile');

subplot(3,2,2);
plot(t,rf);
xlabel('Time (s)');
ylabel('RF (G)');
grid on;
title('RF vs Time');

subplot(3,2,4);
plot(t,grad);
xlabel('Time (s)');
ylabel('Gradient (G/cm)');
grid on;
title('Gradient vs Time');

drawnow;


figure;

subplot(1,3,1);
plot(x,m(1,:));
xlabel('Position (mm)');
ylabel('M_x');
grid on;
title('M_x Profile');

subplot(1,3,2);
plot(x,m(2,:));
xlabel('Position (mm)');
ylabel('M_y');
grid on;
title('M_y Profile');

subplot(1,3,3);
plot(x,m(3,:));
xlabel('Position (mm)');
ylabel('M_z');
grid on;
title('M_z Profile');

drawnow;


% % Plot RF waveform
% f = figure('color','w'); %,'position',[]); 
% hold on;
% plot(t_rf, rf, 'LineWidth', 2);
% plot(t_rf, zeros(size(t_rf)),'-k');
% 
% % Plot RF waveform
% f = figure('color','w'); %,'position',[]); 
% hold on;
% plot(t, rf, 'LineWidth', 2);
% plot(t, zeros(size(t_rf)),'-k');
% yyaxis right; plot(t, grad_ss);

