
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
plot(z,abs(msig));
xlabel('Position (mm)');
ylabel('Signal Magnitude');
grid on;
title('Magnitude Profile');

subplot(3,2,3);
plot(z,angle(msig));
xlabel('Position (mm)');
ylabel('Signal Phase (rad)');
grid on;
title('Phase Profile');

subplot(3,2,5);
plot(z,m(3,:));
xlabel('Position (mm)');
ylabel('Residual M_z');
grid on;
title('M_z Profile');

subplot(3,2,2);
plot(t,rf90);
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
plot(z,m(1,:));
xlabel('Position (mm)');
ylabel('M_x');
grid on;
title('M_x Profile');

subplot(1,3,2);
plot(z,m(2,:));
xlabel('Position (mm)');
ylabel('M_y');
grid on;
title('M_y Profile');

subplot(1,3,3);
plot(z,m(3,:));
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