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
grad_re = [grad_full -0.5*grad_full];
t_re    = [t_full t_full+t_full(end)+dt];


rf90_re = [rf90_full zeros(size(t_full))];



%%% 30 deg flip
% flip angle
alpha = 30; % degrees
% -rf shape- sinc with # of zero-croxxings == TBW
% -rf amplitude- make such that integral of rf waveform gives desired flip
%                   angle (here we want 90 deg filp)
t_rf = 1e3*(t - (tau_rf/2));
rf = sinc(TBW/2 * t_rf);  % sinc(t) = sin(pi*t)/(pi*t)
rf_amp = 1e3 * alpha/(360*gyro*sum(rf));
rf30 = rf_amp * rf;
rf30_full = cat(2,0*grad_rise_up,rf30,0*grad_rise_down);
rf30_re = [rf30_full zeros(size(t_full))];
rf_orig30 = rf;

%%% 10 deg flip
% flip angle
alpha = 10; % degrees
% -rf shape- sinc with # of zero-croxxings == TBW
% -rf amplitude- make such that integral of rf waveform gives desired flip
%                   angle (here we want 90 deg filp)
t_rf = 1e3*(t - (tau_rf/2));
rf = sinc(TBW/2 * t_rf);  % sinc(t) = sin(pi*t)/(pi*t)
rf_amp = 1e3 * alpha/(360*gyro*sum(rf));
rf10 = rf_amp * rf;
rf10_full = cat(2,0*grad_rise_up,rf10,0*grad_rise_down);
rf10_re = [rf10_full zeros(size(t_full))];
rf_orig10 = rf;




%%%% PLOT AND COMPARE FLIP ANGLES
df = 0;

%%% bloch simulation
[msig90,m90] = sliceprofile(rf90_re, grad_re, t_re, ...
    T1, T2, z, df);

%%% bloch simulation
[msig30,m30] = sliceprofile(rf30_re, grad_re, t_re, ...
    T1, T2, z, df);

%%% bloch simulation
[msig10,m10] = sliceprofile(rf10_re, grad_re, t_re, ...
    T1, T2, z, df);


f = figure('color','w','position',[476 498 448 368]);
hold on;
plot(z,abs(msig90),'LineWidth',2);
plot(z,abs(msig30),'LineWidth',2);
plot(z,abs(msig10),'LineWidth',2);
title('Magnitude Slice Profile (df=0Hz)');
legend({'Flip angle = 90deg','Flip angle = 30deg','Flip angle = 10deg'});
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
set(gca,'FontSize',13);

fout = 'plots/p3-part3-flip-angle-comparison.png';
print(f,fout,'-dpng');


%% Slice profile from direct FT of RF pulse waveform

msig_FT_90 = abs(fftshift(fft(rf90)));
msig_FT_10 = abs(fftshift(fft(rf10)));

%% Compute FWHM of slice profiles
hm90 = max(abs(msig90))/2;
fwhmInd90 = find(abs(hm90-abs(msig90))==min(abs(hm90-abs(msig90))));
fwhm90 = z(fwhmInd90(2))-z(fwhmInd90(1));

hm10 = max(abs(msig10))/2;
fwhmInd10 = find(abs(hm10-abs(msig10))==min(abs(hm10-abs(msig10))));
fwhm10 = z(fwhmInd10(2))-z(fwhmInd10(1));

hm90ft = max(abs(msig_FT_90))/2;
fwhmInd90ft = find(abs(hm90ft-abs(msig_FT_90))==min(abs(hm90ft-abs(msig_FT_90))));
fwhm90ft = z(fwhmInd90ft(2))-z(fwhmInd90ft(1));

hm10ft = max(abs(msig_FT_10))/2;
fwhmInd10ft = find(abs(hm10ft-abs(msig_FT_10))==min(abs(hm10ft-abs(msig_FT_10))));
fwhm10ft = z(fwhmInd10ft(2))-z(fwhmInd10ft(1));


%% Plots

f = figure('color','w','position',[476 498 448 368]);
hold on;
plot(z,abs(msig90),'LineWidth',2,'Color','b');
plot(z,abs(msig10),'LineWidth',2,'Color','m');
% yyaxis right; hold on;
plot(z,msig_FT_90,'LineWidth',2,'Color','b','LineStyle',':');
plot(z,msig_FT_10,'LineWidth',2,'Color','m','LineStyle',':');
title('Magnitude Slice Profile (df=0Hz)');
legend({'Bloch, flip=90deg','Bloch, flip=10deg','FT(rf), flip=90deg','FT(rf), flip=10deg'});
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
set(gca,'FontSize',15);

fout = 'plots/p3-part3-blochsim-vs-FT-profile-comparison.png';
print(f,fout,'-dpng');


f = figure('color','w','position',[476 498 448 368]);
hold on;
plot(z(151:251),abs(msig90(151:251)),'LineWidth',2,'Color','b');
plot(z(151:251),abs(msig10(151:251)),'LineWidth',2,'Color','m');
% yyaxis right; hold on;
plot(z(151:251),msig_FT_90(151:251),'LineWidth',2,'Color','b','LineStyle',':');
plot(z(151:251),msig_FT_10(151:251),'LineWidth',2,'Color','m','LineStyle',':');

% plot(z(fwhmInd90),[hm90 hm90],'Color','k');
% plot(z(fwhmInd10),[hm10 hm10],'Color','k');
% plot(z(fwhmInd90ft),[hm90ft hm90ft],'Color','k');
% plot(z(fwhmInd10ft),[hm10ft hm10ft],'Color','k');

title('Magnitude Slice Profile (df=0Hz)');
legend({'Bloch, flip=90deg','Bloch, flip=10deg','FT(rf), flip=90deg','FT(rf), flip=10deg'});
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
set(gca,'FontSize',15);

fout = 'plots/p3-part3-blochsim-vs-FT-profile-comparison-zoom.png';
print(f,fout,'-dpng');



f = figure('color','w','position',[476 224 537 642]);
subplot(211);
hold on;
plot(z,abs(msig90),'LineWidth',2,'Color','b');
% yyaxis right; hold on;
plot(z,msig_FT_90,'LineWidth',2,'Color','b','LineStyle',':');
plot(z(fwhmInd90),[hm90 hm90],'Color','k','LineWidth',2);
plot(z(fwhmInd90ft),[hm90ft hm90ft],'Color','k','LineWidth',2);

title('Magnitude Slice Profile | flip=90deg (df=0Hz)');
legend({'Bloch, flip=90deg','FT(rf), flip=90deg',...
    ['Bloch FWHM=' sprintf('%.2f',fwhm90)],['FT(rf) FWHM=' sprintf('%.2f',fwhm90ft)]});
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
set(gca,'FontSize',13);

subplot(212);
hold on;
plot(z,abs(msig10),'LineWidth',2,'Color','m');
plot(z,msig_FT_10,'LineWidth',2,'Color','m','LineStyle',':');
plot(z(fwhmInd10),[hm10 hm10],'Color','k','LineWidth',2);
plot(z(fwhmInd10ft),[hm10ft hm10ft],'Color','k','LineWidth',2);

legend({'Bloch, flip=10deg','FT(rf), flip=10deg',...
    ['Bloch FWHM=' sprintf('%.2f',fwhm10)],['FT(rf) FWHM=' sprintf('%.2f',fwhm10ft)]});
title('Magnitude Slice Profile | flip=10deg (df=0Hz)');
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
set(gca,'FontSize',13);

fout = 'plots/p3-part3-blochsim-vs-FT-profile-comparison_subplots+fwhm.png';
print(f,fout,'-dpng');





%% Change T2 to 2ms

%%%% PLOT AND COMPARE FLIP ANGLES
df = 0;

%%% bloch simulation
[msig90_t2_2,m90_t2_2] = sliceprofile(rf90_re, grad_re, t_re, ...
    T1, 2, z, df);

%%% bloch simulation
[msig30_t2_2,m30_t2_2] = sliceprofile(rf30_re, grad_re, t_re, ...
    T1, 2, z, df);

%%% bloch simulation
[msig10_t2_2,m10_t2_2] = sliceprofile(rf10_re, grad_re, t_re, ...
    T1, 2, z, df);



f = figure('color','w','position',[476 464 580 402]);
hold on;
p1=plot(z,abs(msig90),'LineWidth',2);
plot(z,abs(msig90_t2_2),'LineWidth',2,'Color',p1.Color,'LineStyle',':');
p2=plot(z,abs(msig10),'LineWidth',2);
plot(z,abs(msig10_t2_2),'LineWidth',2,'Color',p2.Color,'LineStyle',':');

title('Magnitude Slice Profile (df=0Hz)');
legend({'Flip=90deg, T2=100ms','Flip=90deg, T2=2ms','Flip=30deg, T2=100ms','Flip=30deg, T2=2ms'});
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
set(gca,'FontSize',15);

fout = 'plots/p3-part3-T2-comparison.png';
print(f,fout,'-dpng');


%% P3 Part 4: 
% For one of the pulses you generated, simulate mx, my, mz and mxy at the 
% point in time between the slice selective gradient and the slice 
% rephasing gradient. Explain the purpose of the slice rephasing gradient.



