function f = make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, df, saveplots )
% function f = make_p3part2_figures_re(z, rf90_re, grad_re, m, msig, TBW, tau_rf, t_re, df, saveplots )
% 
% Wrapper function to generate plots/figures for problem 3 part 2 
% (See p2part2.m)
% RJ | 10-03-2023

if nargin<10
    saveplots=0;
end

% Plot RF waveform
f = figure('color','w','position',[476 96 562 770]); 

subplot(211);
hold on;
plot(t_re*1e3, rf90_re*1e2, 'LineWidth', 2);
plot(t_re*1e3, zeros(size(t_re)),'-k');
set(gca,'FontSize',13);
xlabel('Time (ms)');
ylabel('RF amplitude (\muT)');
title({'P3.2 | RF pulse waveform', ...
    ['Off-resonance freq = ' sprintf('%g',df) ' Hz']});

subplot(212);
hold on;
plot(t_re*1e3, grad_re*1e-1, 'LineWidth', 2);
plot(t_re*1e3, zeros(size(t_re)),'-k');
set(gca,'FontSize',13);
xlabel('Time (ms)');
ylabel('G_{ss} amplitude (mT)');
title({'P3.2 | G_{ss} waveform'}); % ...
%     ['TBW=' num2str(TBW) ', \tau_{RF}=' num2str(tau_rf)]});

if saveplots
fout = ['plots/p3-part2-rf+grad_df' sprintf('%g',df) 'hz.png'];
print(f,fout,'-dpng');
end

f = figure('color','w','position',[470 51 493 801]);
subplot(311);
hold on;
plot(z, m(1,:), 'LineWidth', 2);
plot(z, m(2,:), 'LineWidth', 2);
plot(z, m(3,:), 'LineWidth', 2);
plot(z, zeros(size(z)),'-k');
ylim([-1 1]);
set(gca,'FontSize',13);
xlabel('Position (mm)');
ylabel('Residual Magnetization');
lgn = legend({'M_x','M_y','M_z'},'Position',[0.707662271805274 0.727282658293597 0.114604462474645 0.0774032459425719]);
title({'P3.2 | M profiles'}); % ...
%     ['TBW=' num2str(TBW) ', \tau_{RF}=' num2str(tau_rf)]});

subplot(312);
hold on;
plot(z, abs(msig), 'LineWidth', 2);
plot(z, zeros(size(z)),'-k');
set(gca,'FontSize',13);
ylim([0 1]);
xlabel('Position (mm)');
ylabel('Residual M_{xy} - Signal Magnitude');
title({'P3.2 | M_{xy} profile (abs)'}); % ...
%     ['TBW=' num2str(TBW) ', \tau_{RF}=' num2str(tau_rf)]});

subplot(313);
hold on;
plot(z, angle(msig), 'LineWidth', 2);
plot(z, zeros(size(z)),'-k');
set(gca,'FontSize',13);
ylim([-pi pi]);
xlabel('Position (mm)');
ylabel('Residual M_{xy} - Signal Phase');
title({'P3.2 | M_{xy} profile (angle)'}); % ...
%     ['TBW=' num2str(TBW) ', \tau_{RF}=' num2str(tau_rf)]});

if saveplots
% fout = 'plots/p3-part2-rf+grad.png';
fout = ['plots/p3-part2-mprofiles_df' sprintf('%g',df) 'hz.png'];
print(f,fout,'-dpng');
end


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
% plot(t_re,rf90);
% xlabel('Time (s)');
% ylabel('RF (G)');
% grid on;
% title('RF vs Time');
% 
% subplot(3,2,4);
% plot(t_re,grad_re);
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