

t0 = [-1:.001:1];
x= [-5:.1:5];

% [msig,m]=sliceprofile(.05*sinc(1000*t-3), 0.1*ones(size(t)), t, 600, 100, x, 0);
%  [msig,m]=sliceprofile(rf,grad,t,T1,T2,pos,df)

t = [-1:.001:1];


TBW = 4;

rf = sinc(TBW*t);
rf_amp = deg2rad(90) / (42.58*10^6 * sum(rf));
rf = rf * rf_amp;

figure; plot(t,rf);

df = 0;
x= [-100:.1:100]*1e-3;
T1=1000;
T2=100;
grad_amp = 18.788 * 1e-3 * 1e4 * 1e-2;

t = [0:.001:2];

[msig,m] = sliceprofile(rf, grad_amp*ones(size(t)), t, T1, T2, x, df);

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
