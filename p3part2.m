% -----------------------------------------
% BME599 F23 | HW1 P3(2) - Slice Profile Simulation
% Bloch Equation Simulation, bSSFP sequence
% -----------------------------------------
% Robert Jones | 09-24-2023
% -----------------------------------------

% t = [0:.0001:.006]; % s
% x= [-20:.1:20];     % mm
% [msig,m]=sliceprofile(.05*sinc(1000*t-3), 0.1*ones(size(t)), t, 600, 100, x, 0);
%  [msig,m]=sliceprofile( rf,  grad,  t,  T1, T2,  pos,  df)
%                  UNITS: 

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
rf_full = cat(2,0*grad_rise_up,rf90,0*grad_rise_down);

% %%% Add rephasing gradient (no grad ramping)
% grad_re = [grad -0.5*grad];
% rf90_re = [rf90 zeros(size(t))];
% t_re    = [t t+t(end)+dt];

%%% Add rephasing gradient (with grad ramping)
grad_re = [grad_full -0.5*grad_full];
rf90_re = [rf_full zeros(size(t_full))];
t_re    = [t_full t_full+t_full(end)+dt];

%%% off-resonance
df = 0; % Hz

%%% bloch simulation
[msig,m] = sliceprofile(rf90_re, grad_re, t_re, ...
    T1, T2, z, df);

%%% make plots
SavePlots=1;
% make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, plotTag );
make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, df, SavePlots )

%%% off-resonance
df = 200; % Hz

%%% bloch simulation
[msig,m] = sliceprofile(rf90_re, grad_re, t_re, ...
    T1, T2, z, df);

%%% make plots
SavePlots=1;
% make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, plotTag );
make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, df, SavePlots );




