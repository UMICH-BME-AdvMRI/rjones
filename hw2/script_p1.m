%%  BME599 - HW2 P1
%   Extended Phase Graphs
%   Robert Jones | 10-24-2023

%%%%  PLEASE NOTE:          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This code uses functions from the Stanford EPG Matlab toolbox. The code
%  can be found here: https://web.stanford.edu/~bah/software/epg/


%% Setup
addpath('epg_funcs');
addpath('blochsim_utils');

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
    200, 50;
    500, 100;
    800, 150;
    1200, 200;
    1500, 300;
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
print(f,[plotdir filesep 'p1a-i.png'],'-dpng');

part1a1.results = results;
part1a1.params = params;


%% Part (a)-ii
% Set parameters for 
params.flip = deg2rad(120);   % flip angle (radians)

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
title({ 'P1a(ii) | FSE EPG simulation', sprintf('etl=%d, esp=%dms, flip=%gdeg', ...
    params.etl, params.esp, rad2deg(params.flip)) });
set(gca,'FontSize',16);
print(f,[plotdir filesep 'p1a-ii.png'],'-dpng');

part1a2.results = results;
part1a2.params = params;


%% Part (a)-iii
% Set parameters for 
params.flip = deg2rad(60);   % flip angle (radians)

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
title({ 'P1a(iii) | FSE EPG simulation', sprintf('etl=%d, esp=%dms, flip=%gdeg', ...
    params.etl, params.esp, rad2deg(params.flip)) });
set(gca,'FontSize',16);
print(f,[plotdir filesep 'p1a-iii.png'],'-dpng');

part1a3.results = results;
part1a3.params = params;


%% Part (b)
%   Plot 4 contour plots signals vs T1 and T2 
%   for the 6th, 16th, 32nd, and 48th echoes using various α

% Set parameters for 
params.etl = 64;    % echo train length
params.esp = 5;     % echo spacing
params.echotimes = (1:params.etl)*params.esp;  %times of echos

results = [];

T1s = 200:100:1500;
T2s = 50:30:300;
echos = [6,16,32,48];
flipangles = 60:60:180; % 30:30:180;

signals = zeros(length(T2s),length(T1s),length(echos),length(flipangles));

t1labels = {}; t2labels = {};
for ff=1:length(flipangles)
    flipang = flipangles(ff);
    disp(flipang);
    for ind2=1:length(T2s)
%         ind2inv = length(T2s)-ind2+1;
        t2 = T2s(ind2);
        t2labels{ind2} = num2str(t2);
        for ind1=1:length(T1s)
            t1 = T1s(ind1);
            t1labels{ind1} = num2str(t1);

            [s1,pd1,p1] = epg_cpmg(deg2rad(flipang), params.etl, t1, t2, params.esp);
        	signals(ind2,ind1,:,ff) = s1(echos);
            
        end
    end
end


%%%%%%%% CONTOUR PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make different plots, one for each flip angle
%   Each with 4 sublots, one for each echo 
for ff=1:length(flipangles)
    flipang = flipangles(ff);
    disp(flipang);
    dat = squeeze(abs(signals(:,:,:,ff)));
    figure1 = createfigure_p1b(T1s, T2s, dat(:,:,1), dat(:,:,2), dat(:,:,3), dat(:,:,4), flipang);
    shg;
    print(figure1,[plotdir filesep 'p1b_flip' num2str(flipang) '_contour.png'],'-dpng');
end



%% Removed code

if 0 == 1

    %%%%%%%% IMAGE PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make 6 different plots, one for each flip angle
    %   Each with 4 sublots, one for each echo 
    for ff=1:length(flipangles)
        flipang = flipangles(ff);
        disp(flipang);
    
        f = figure('position',[60 586 1425 267],'color','w');
        for ee=1:length(echos)
            subplot(1,4,ee);
            dat = squeeze(abs(signals(:,:,ee,ff)));
            dat = flipud(dat);
            bigdat = imresize(dat,100,'nearest');
            imagesc((dat));
    %         surf(dat);
    %         xlabel('T1 (ms)');
    %         ylabel('T2 (ms)');
    %         set(gca,'XTick',1:2:length(T1s));
    %         set(gca,'XTickLabel',t1labels(1:2:end));
    % %         set(gca,'YTick',1:1:length(T2s));
    % %         set(gca,'YTickLabel',t2labels);
    %         set(gca,'YTick',1:1:length(T2s));
    %         set(gca,'YTickLabel',t2labels(end:-1:1));
            title(sprintf('Echo %d, Flip = %ddeg',echos(ee),flipang));
    %         colorbar; 
    %         colormap gray;
    %         caxis([0 1]);
            set(gca,'FontSize',12);
            drawnow;
        end
    %     print(f,[plotdir filesep 'new_p1b_flip' num2str(flipang) '.png'],'-dpng');
    end
    
    %%%%%%%% CONTOUR PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make different plots, one for each flip angle
    %   Each with 4 sublots, one for each echo 
    for ff=1:length(flipangles)
        flipang = flipangles(ff);
        disp(flipang);
    
        f = figure('position',[60 586 1425 267],'color','w');
        for ee=1:length(echos)
            subplot(1,4,ee);
    %         subplot(2,4,ee);
            ax=gca;
            dat = squeeze(abs(signals(:,:,ee,ff)));
            [c,h] = contour(T1s,T2s,dat);
    %         set(gca, 'XDir','reverse')
            h.LineWidth = 2;
            xlabel('T1 (ms)');
            ylabel('T2 (ms)');
            title(sprintf('Echo %d, Flip = %ddeg',echos(ee),flipang));
            set(gca,'FontSize',12);
    
    %         subplot(2,4,ee+4);
    %         imagesc(imresize(dat,100));set(gca,'colormap',jet);
    %         drawnow;
        end
    %     print(f,[plotdir filesep 'new_p1b_flip' num2str(flipang) '_contour.png'],'-dpng');
    end
    
end


if 0 == 1
    % Make 4 different plots, one for each echo
    %   Each with 6 sublots, one for each flip 
    t1labels = {};
    for ii=1:length(T1s)
        t1labels{ii} = sprintf('%d',T1s(ii));
    end
    t2labels = {};
    for ii=1:length(T2s)
        t2labels{ii} = sprintf('%d',T2s(ii));
    end
    
    for ee=1:length(echos)
        f = figure('position',[60 586 1425 267],'color','w');
        for ff=1:length(flipangles)
            flipang = flipangles(ff);
            disp(flipang);
            subplot(1,6,ff);
            dat = squeeze(abs(signals(:,:,ee,ff)));
            imagesc(flipud(dat.'));
            xlabel('T1 (ms)');
            ylabel('T2 (ms)');
            set(gca,'XTick',1:2:length(T1s));
            set(gca,'XTickLabel',t1labels(1:2:end));
    %         set(gca,'YTick',1:1:length(T2s));
    %         set(gca,'YTickLabel',t2labels);
            set(gca,'YTick',1:1:length(T2s));
            set(gca,'YTickLabel',t2labels(end:-1:1));
            title(sprintf('Echo %d, Flip = %ddeg',echos(ee),flipang));
    %         colorbar; 
            colormap gray;
            caxis([0 1]);
            set(gca,'FontSize',12);
            drawnow;
        end
        print(f,[plotdir filesep 'p1b_echo' num2str(echos(ee)) '.png'],'-dpng');
    end
end


% Make one big plot showing signal images
if 0 == 1
    f = figure('position',[60 37 1425 816],'color','w');
    
    for ff=1:length(flipangles)
        flipang = flipangles(ff);
        disp(flipang);
    
        
        for ee=1:length(echos)
            ppind=4*(ff-1)+ee;
            subplot(6,4,ppind);
            dat = squeeze(abs(signals(:,:,ee,ff)));
            imagesc(flipud(dat.'));
            xlabel('T1 (ms)');
            ylabel('T2 (ms)');
            set(gca,'XTick',1:2:length(T1s));
            set(gca,'XTickLabel',t1labels(1:2:end));
    %         set(gca,'YTick',1:1:length(T2s));
    %         set(gca,'YTickLabel',t2labels);
            set(gca,'YTick',1:1:length(T2s));
            set(gca,'YTickLabel',t2labels(end:-1:1));
            title(sprintf('Echo %d, Flip = %ddeg',echos(ee),flipang));
    %         colorbar; 
            colormap gray;
            caxis([0 1]);
            set(gca,'FontSize',12);
            drawnow;
        end
    end

end

% % Run EPG simulations
% results.s = [];
% results.phasediag = [];
% results.P = [];
% results.t1t2 = [];
% results.labels = {};
% for tind=1:size(params.T1T2s,1)
%     currt = params.T1T2s(tind,:);
%     disp(currt);
%     currt1=currt(1); currt2=currt(2);
%     [s1,pd1,p1] = epg_cpmg(params.flip, params.etl, currt1, currt2, params.esp);	% IDEAL, CPMG
% 
%     results.t1t2(tind,:) = currt;
%     results.s = cat(2,results.s,reshape(s1,[],1));
%     results.phasediag = cat(3,results.phasediag,pd1);
%     results.P = cat(3,results.P,p1);
%     results.labels{tind} = sprintf('T1=%dms, T2=%dms',currt(1),currt(2));
% end




% % % % FROM TOOLBOX: epg_echotrain.m
%
% s2 = epg_cpmg(pi*i,30,1000,200,10);	% IDEAL, non CPMG
% s3 = epg_cpmg(2/3*pi,30,1000,200,10);	% 120 CPMG
% flips = 2/3*pi*ones(1,30);		% Constant flip angles
% s4 = epg_cpmg(flips,30,1000,200,10);	% 120 without 1st pulse correction.
% s5 = epg_cpmg(2/3*pi*i,30,1000,200,10);	% 120 non CPMG
% t = [1:30]*10;	% Echo times.
% plot(t,abs([s1(:) s2(:) s3(:) s4(:) s5(:)]));
% legend('180 CPMG','180 Non-CPMG','120 CPMG','120 CPMG no corr','120 Non-CPMG');

