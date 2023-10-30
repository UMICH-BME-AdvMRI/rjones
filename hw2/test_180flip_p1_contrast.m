%% Setup
addpath('epg_funcs');

clear
% close all;
% clc

% To save 
plotdir = '/Users/robertjones/Desktop/F23/599/hw/hw2/code2upload/plots_p1';
if ~exist(plotdir,'dir'), mkdir(plotdir); end


%% Part (a)-i
% Set parameters for 
params.etl = 64;    % echo train length
params.esp = 5;     % echo spacing
params.echotimes = (1:params.etl)*params.esp;  %times of echos
params.T1T2s = [
    200, 100;
    500, 100;
    800, 100;
    1200, 100;
    1500, 100;
    ];              % 5 different T1, T2 values to test
params.flip = pi;   % flip angle (radians)

% Run EPG simulations
results.s = [];
results.phasediag = [];
results.P = [];
results.t1t2 = [];
results.labels = {};
for tind=1:size(params.T1T2s,1)
    currt = params.T1T2s(tind,:);
    disp(currt);
    currt1=currt(1); currt2=currt(2);
    [s1,pd1,p1] = epg_cpmg(params.flip, params.etl, currt1, currt2, params.esp);	% IDEAL, CPMG

    results.t1t2(tind,:) = currt;
    results.s = cat(2,results.s,reshape(s1,[],1));
    results.phasediag = cat(3,results.phasediag,pd1);
    results.P = cat(3,results.P,p1);
    results.labels{tind} = sprintf('T1=%dms, T2=%dms',currt(1),currt(2));
end

% Plot the echo amplitudes
f = figure('position',[832 344 673 508],'color','w'); 
hold on;
for tind=1:size(params.T1T2s,1)
    currt = params.T1T2s(tind,:);
    y = results.s(:,tind);
    plot(params.echotimes,abs(y),'LineWidth',2);
    drawnow;
%     pause(1);
end
lgn = legend(results.labels,'location','best');
xlabel('Echo time (ms)');
ylabel('Echo amplitude');
xlim([params.echotimes(1) params.echotimes(end)]);
title({ 'P1a(i) | FSE EPG simulation', sprintf('etl=%d, esp=%dms, flip=%gdeg', ...
    params.etl, params.esp, rad2deg(params.flip)) });
set(gca,'FontSize',16);
print(f,[plotdir filesep 'varying-T1_constant-T2_flip180deg.png'],'-dpng');





%% 

params.etl = 64;    % echo train length
params.esp = 5;     % echo spacing
params.echotimes = (1:params.etl)*params.esp;  %times of echos
params.T1T2s = [
    1000, 10;
    1000, 40;
    1000, 70;
    1000, 100;
    1000, 130;
    ];              % 5 different T1, T2 values to test
params.flip = pi;   % flip angle (radians)

% Run EPG simulations
results.s = [];
results.phasediag = [];
results.P = [];
results.t1t2 = [];
results.labels = {};
for tind=1:size(params.T1T2s,1)
    currt = params.T1T2s(tind,:);
    disp(currt);
    currt1=currt(1); currt2=currt(2);
    [s1,pd1,p1] = epg_cpmg(params.flip, params.etl, currt1, currt2, params.esp);	% IDEAL, CPMG

    results.t1t2(tind,:) = currt;
    results.s = cat(2,results.s,reshape(s1,[],1));
    results.phasediag = cat(3,results.phasediag,pd1);
    results.P = cat(3,results.P,p1);
    results.labels{tind} = sprintf('T1=%dms, T2=%dms',currt(1),currt(2));
end

% Plot the echo amplitudes
f = figure('position',[832 344 673 508],'color','w'); 
hold on;
for tind=1:size(params.T1T2s,1)
    currt = params.T1T2s(tind,:);
    y = results.s(:,tind);
    plot(params.echotimes,abs(y),'LineWidth',2);
    drawnow;
%     pause(1);
end
lgn = legend(results.labels,'location','best');
xlabel('Echo time (ms)');
ylabel('Echo amplitude');
xlim([params.echotimes(1) params.echotimes(end)]);
title({ 'P1a(i) | FSE EPG simulation', sprintf('etl=%d, esp=%dms, flip=%gdeg', ...
    params.etl, params.esp, rad2deg(params.flip)) });
set(gca,'FontSize',16);
print(f,[plotdir filesep 'constant-T2_varying-T1_flip180deg.png'],'-dpng');

