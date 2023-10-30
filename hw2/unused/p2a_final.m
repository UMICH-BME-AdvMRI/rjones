%%  BME599 - HW2 P2
%   Single and Multiple Spin Echo Sequences
%   Robert Jones | 10-24-2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  PLEASE NOTE:          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code uses functions from the Stanford EPG and BlochSim Matlab toolboxes. 
% The EPG code can be found here: https://web.stanford.edu/~bah/software/epg/
% The BlochSim code can be found here: http://mrsrl.stanford.edu/~brian/bloch/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
addpath('epg_funcs');
addpath('/Users/robertjones/Desktop/F23/599/hw/hw1/blochsim_matlab/utils');

clear
% close all;
% clc

% To save 
plotdir = '/Users/robertjones/Desktop/F23/599/hw/hw2/plots_p2';
if ~exist(plotdir,'dir'), mkdir(plotdir); end

load('brain_maps.mat');
dims = size(T1map);

f = figure('position',[175 509 1051 265],'color','w');
subplot(131); 
imshow(M0map,[]); title('M0 map'); set(gca,'FontSize',16);
subplot(132);
imshow(T1map,[]); title('T1 map'); set(gca,'FontSize',16);
subplot(133);
imshow(T2map,[]); title('T2 map'); set(gca,'FontSize',16);
print(f,[plotdir filesep 'input_brain_maps.png'],'-dpng');


% Make binary mask of voxels to process
mask0 = M0map==0; mask1 = T1map==0; mask2 = T2map==0;
mask3 = mask1 | mask2 | mask0;
mask = ~mask3;

%% Part (a) Single spin echo
% Simulate the following contrasts:
%    T1-weighted, T2-weighted, and proton density weighted

warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T1w - short TR short TE - TR/TE = 500/15 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TE = 15;
TR = 500;
dfreq = 0;

t1w = zeros(dims);
for xx=1:dims(1)
    if mod(xx,25)==0, disp(xx); end
    for yy=1:dims(2)
        if mask(xx,yy)==1
            m0val = M0map(xx,yy);
            t1val = T1map(xx,yy);
            t2val = T2map(xx,yy);
    
            [Msig,Mss] = sesignal(t1val,t2val,TE,TR,dfreq,m0val);
            t1w(xx,yy) = abs(Msig);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T2w - long TR long TE - TR/TE = 6000/100 ms %4000/90 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TE = 100;
TR = 6000;
dfreq = 0;

t2w = zeros(dims);
for xx=1:dims(1)
    if mod(xx,25)==0, disp(xx); end
    for yy=1:dims(2)
        if mask(xx,yy)==1
            m0val = M0map(xx,yy);
            t1val = T1map(xx,yy);
            t2val = T2map(xx,yy);
    
            [Msig,Mss] = sesignal(t1val,t2val,TE,TR,dfreq,m0val);
            t2w(xx,yy) = abs(Msig);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PD - long TR short TE - TR/TE = 6000/15 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TE = 15;
TR = 6000;
dfreq = 0;

pd = zeros(dims);
for xx=1:dims(1)
    if mod(xx,25)==0, disp(xx); end
    for yy=1:dims(2)
        if mask(xx,yy)==1
            m0val = M0map(xx,yy);
            t1val = T1map(xx,yy);
            t2val = T2map(xx,yy);
    
            [Msig,Mss] = sesignal(t1val,t2val,TE,TR,dfreq,m0val);
            pd(xx,yy) = abs(Msig);
        end
    end
end

% Display t1w, t2w, pd simulated contrasts
f = figure('color','w','position',[312 494 1040 344]);
subplot(131);
imshow(t1w,[]);
title({'T1w',['TR/TE = 500/15 ms']});
set(gca,'FontSize',16);
subplot(132);
imshow(t2w,[]);
title({'T2w',['TR/TE = 6000/100 ms']});
set(gca,'FontSize',16);
subplot(133);
imshow(pd,[]);
% title({'PD',['TR/TE = 4000/15 ms']});
title({'PD',['TR/TE = 6000/15 ms']});
set(gca,'FontSize',16);
print(f,[plotdir filesep 'p2a_t1w,t2w,pd.png'],'-dpng');



%% Part b(i) - using EPG cpmg modified for multiple TRs
TR = 3000;
esp = 5;
etl = 32;
numTRs = 5;
flipangle = pi;
m0val = 1;
echotimes = esp*(1:etl);

T1s = [1000,1000,2000,2000];
T2s = [50,100,50,100];

ntrials = length(T1s);

Msss = zeros(etl,3,ntrials);
Msigs = zeros(etl,ntrials);

timelabels = cell(ntrials,1);
for nn=1:ntrials
    t1val = T1s(nn);
    t2val = T2s(nn);
    [s,phasediag,P] = epg_cpmg_multi_tr(flipangle,etl,t1val,t2val,esp,TR,numTRs,m0val,0);
    Msigs(:,nn) = abs(s);
    timelabels{nn} = sprintf('TE/TR=%d/%d ms',t1val,t2val);
end

f = figure('color','w'); %,'position',[]); 
plot(echotimes,Msigs,'LineWidth',2);
legend(timelabels);
xlabel('Echo times (ms)');
ylabel('Signal (Transverse magnetization)');
set(gca,'FontSize',15);
title('P2b(i) | FSE signal after 5 TRs');
print(f,[plotdir filesep 'p2b-i.png'],'-dpng');

%%%%% Part b(i) - using bloch sim FSE
dfreq = 0;
Msss_bloch = zeros(etl,3,ntrials);
Msigs_bloch = zeros(etl,ntrials);
for nn=1:ntrials
    t1val = T1s(nn);
    t2val = T2s(nn);  
    [Msig,Mss] = fsesignal_multi(t1val,t2val,esp,TR,dfreq,etl,m0val,numTRs);
    Msigs_bloch(:,nn) = abs(Msig);
    Msss_bloch(:,:,nn) = Mss.';
end
f = figure('color','w'); %,'position',[]); 
hold on;
plot(echotimes,Msigs,'LineWidth',2);
plot(echotimes,Msigs_bloch,'LineWidth',2,'LineStyle','--');
legend('TR/TE = 1000/50 ms','TR/TE = 1000/100 ms','TR/TE = 2000/50 ms',...
    'TR/TE = 1000/50 ms (bloch)','TR/TE = 1000/100 ms (bloch)','TR/TE = 2000/50 ms (bloch)');
xlabel('Echo times (ms)');
ylabel('Signal (Transverse magnetization)');
set(gca,'FontSize',15);
title('P2b(i) | FSE signal after 5 TRs');
print(f,[plotdir filesep 'p2b-i_bloch-v-epg-comparison.png'],'-dpng');



%% Part b(ii)

TR = 3000;
esp = 5;
etl = 32;
echotimes = esp*(1:etl);
numTRs = 1;
flipangle = pi;
brain_mese = zeros([dims,etl]);
brain_mese_epg = zeros([dims,etl]);
tic;
for xx=1:dims(1)
    if mod(xx,25)==0, disp(xx); end
    for yy=1:dims(2)
        if mask(xx,yy)==1
            t1val = T1map(xx,yy);
            t2val = T2map(xx,yy);
            m0val = M0map(xx,yy);
            % Use epg simulation
            [s,~,~] = epg_cpmg_multi_tr(flipangle,etl,t1val,t2val,esp,TR,numTRs,m0val,0);
            brain_mese_epg(xx,yy,:) = (s);
            % Use blochsim simulation
            [Msig,Mss] = fsesignal_multi(t1val,t2val,esp,TR,dfreq,etl,m0val,numTRs);
            brain_mese(xx,yy,:) = (Msig);
        end
    end
end
toc

%%%% Fill k-space -  TE_eff = 80 ms
%    - Use the middle (16th) echo to fill kspace center line
%    - So, just fill k-space from top->bottom with echos 1->etl

% Desired TE_eff
te_eff = 80; %ms
% number of kspace lines to use per echo
klines_per_echo = dims(1)/etl;

% Index of echo to fill k-space center
central_echo_ind = echotimes==te_eff;
% the central line of k-space
dc_offset = dims(1)/2+1;

% Visualize k-space filling order:
k_filling_order = 1:etl;
k_fill_image = zeros(dims);
figure('position',[1005 251 484 357]); 
hold on; 
for ee=1:etl
    kstart = (etl-ee)*klines_per_echo+1;
    kend   = (etl-ee+1)*klines_per_echo;
    k_fill_image(kstart:kend,:) = ee;
    imshow(k_fill_image,[0 etl]);
    title(sprintf('echo %d, ky=%d,%d',ee,kstart,kend));
    colormap jet; 
    set(gca,'FontSize',16);
    set(gcf,'position',[1077 230 412 378]);
%     shg; pause;
end
title(['TE_{eff} = ' num2str(te_eff) 'ms, etl = ' num2str(etl)]);
cbar=colorbar;
cbar.Title.String = 'Echo #';
xlabel('kx'); ylabel('ky');
k_fill_image_80 = k_fill_image;

% Generate k-space
brain_kspace = zeros(dims);
for ee=1:etl
    kstart = (etl-ee)*klines_per_echo+1;
    kend   = (etl-ee+1)*klines_per_echo;
    freqim = fftshift(fft2(brain_mese(:,:,ee)));
    brain_kspace(kstart:kend,:) = freqim(kstart:kend,:);
end
brain_sim = abs(ifft2(fftshift(brain_kspace)));
brain_te80_etl32 = brain_sim;
figure; imshow(brain_te80_etl32,[]);



%%%% Fill k-space -  TE_eff = 40 ms
%    - Use the 8th echo to fill kspace center line

% Desired TE_eff
te_eff = 40; %ms
% Index of echo to fill k-space center
central_echo_ind = find(echotimes==te_eff);

% Visualize k-space filling order:
k_filling_order = 1:etl;
k_fill_image = zeros(dims);
figure('position',[1005 251 484 357]); 
hold on; 
for ee=1:etl
    ee_offset = ee+central_echo_ind;
    kstart = (etl-ee_offset)*klines_per_echo+1;
    kend   = (etl-ee_offset+1)*klines_per_echo;
    if kstart<1, kstart=kstart+dims(1); end
    if kend<1, kend=kend+dims(1); end
    k_fill_image(kstart:kend,:) = ee;
    imshow(k_fill_image,[0 etl]);
    title(sprintf('echo %d, ky=%d,%d',ee,kstart,kend));
    colormap jet; 
    set(gca,'FontSize',16);
    set(gcf,'position',[1077 230 412 378]);
%     shg; pause;
end
title(['TE_{eff} = ' num2str(te_eff) 'ms, etl = ' num2str(etl)]);
cbar=colorbar;
cbar.Title.String = 'Echo #';
xlabel('kx'); ylabel('ky');
k_fill_image_40 = k_fill_image;

% Generate k-space
brain_kspace = zeros(dims);
for ee=1:etl
    ee_offset = ee+central_echo_ind;
    kstart = (etl-ee_offset)*klines_per_echo+1;
    kend   = (etl-ee_offset+1)*klines_per_echo;
    if kstart<1, kstart=kstart+dims(1); end
    if kend<1, kend=kend+dims(1); end
    freqim = fftshift(fft2(brain_mese(:,:,ee)));
    brain_kspace(kstart:kend,:) = freqim(kstart:kend,:);
end
brain_sim = abs(ifft2(fftshift(brain_kspace)));
brain_te40_etl32 = brain_sim;


%%%% Fill k-space -  TE_eff = 120 ms
%    - Use the 24th echo to fill kspace center line

% Desired TE_eff
te_eff = 120; %ms
% Index of echo to fill k-space center
central_echo_ind = find(echotimes==te_eff);

% Visualize k-space filling order:
k_filling_order = 1:etl;
k_fill_image = zeros(dims);
figure('position',[1005 251 484 357]); 
hold on; 
for ee=1:etl
    ee_offset = ee+central_echo_ind;
    kstart = (etl-ee_offset)*klines_per_echo+1;
    kend   = (etl-ee_offset+1)*klines_per_echo;
    if kstart<1, kstart=kstart+dims(1); end
    if kend<1, kend=kend+dims(1); end
    k_fill_image(kstart:kend,:) = ee;
    imshow(k_fill_image,[0 etl]);
    title(sprintf('echo %d, ky=%d,%d',ee,kstart,kend));
    colormap jet; 
    set(gca,'FontSize',16);
    set(gcf,'position',[1077 230 412 378]);
%     shg; pause;
end
title(['TE_{eff} = ' num2str(te_eff) 'ms, etl = ' num2str(etl)]);
cbar=colorbar;
cbar.Title.String = 'Echo #';
xlabel('kx'); ylabel('ky');
k_fill_image_120 = k_fill_image;

% Generate k-space
brain_kspace = zeros(dims);
for ee=1:etl
    ee_offset = ee+central_echo_ind;
    kstart = (etl-ee_offset)*klines_per_echo+1;
    kend   = (etl-ee_offset+1)*klines_per_echo;
    if kstart<1, kstart=kstart+dims(1); end
    if kend<1, kend=kend+dims(1); end
    freqim = fftshift(fft2(brain_mese(:,:,ee)));
    brain_kspace(kstart:kend,:) = freqim(kstart:kend,:);
end
brain_sim = abs(ifft2(fftshift(brain_kspace)));
brain_te120_etl32 = brain_sim;




%%% PLOT IMAGES
f = figure('color','w','position',[476 572 931 294]);
subplot(131);
imshow(brain_te40_etl32,[]); %[0 1]);
title('FSE sim | TE_{eff} = 40 ms');
set(gca,'FontSize',16);
subplot(132);
imshow(brain_te80_etl32,[]); %[0 1]);
title('FSE sim | TE_{eff} = 80 ms');
set(gca,'FontSize',16);
subplot(133);
imshow(brain_te120_etl32,[]); %[0 1]);
title('FSE sim | TE_{eff} = 120 ms');
set(gca,'FontSize',16);
print(f,[plotdir filesep 'p2b-ii+iii.png'],'-dpng');



%% Part b(iv)
% Set parameters to use
TR = 3000;
esp = 5;
numTRs = 1;
flipangle = pi;
etls = [16,32,64,128];
te_eff = 80;
% echoind = find(echotimes==te_eff);

% To store MESE images for different ETLs
brain_mese_etls = cell(length(etls),1);

% Store results from ETL==32, and skip in for loop below
brain_mese_etls{find(etls==32),1} = brain_mese;

for etlind = 1:length(etls)
    etl = etls(etlind);
    if etl~=32
        fprintf(' -on etl = %d\n',etl);
        brain_mese = zeros([dims,etl]);
        for xx=1:dims(1)
            if mod(xx,10)==0, fprintf(' x=%d\n',xx); end
            for yy=1:dims(2)
                if mask(xx,yy)==1
                    t1val = T1map(xx,yy);
                    t2val = T2map(xx,yy);
                    m0val = M0map(xx,yy);
                    [Msig,Mss] = fsesignal_multi(t1val,t2val,esp,TR,dfreq,etl,m0val,numTRs);
                    brain_mese(xx,yy,:) = (Msig);
                end
            end
        end
        brain_mese_etls{find(etls==etl),1} = brain_mese;
    end
end

% To store images with TE_eff=80ms and different ETLs
brain_sim_te80 = zeros([dims,length(etl)]);

% Simulate images with TE_eff=80ms and different ETLs
for tt=1:length(etls)
    etl = etls(tt);
    echotimes = (1:etl)*esp;
    central_echo_ind = find(echotimes==te_eff);
    klines_per_echo = dims(1)/etl;
    
    % Set kspace PE line filling order based on ETL
    switch etl
        case 16
            klineorder = [8:16,1:7];
        case 32
            klineorder = etl:-1:1;
        case 64
            klineorder = [48:etl,1:15,16:47];
        case 128
            klineorder = [80:128,1:15,16:79];
    end

    % Visualize k-space filling order:
    k_fill_image = zeros(dims);
    k_echo_lines = zeros(etl,2);
    figure('position',[1005 251 484 357]); 
    hold on; 
    for ee=1:etl
        kstart = (ee-1)*klines_per_echo+1;
        kend   = kstart+klines_per_echo-1;
        k_fill_image(kstart:kend,:) = klineorder(ee);
        k_echo_lines(klineorder(ee),:) = [kstart,kend];
        imshow(k_fill_image,[0 etl]);
        title(sprintf('echo %d, ky=%d,%d',klineorder(ee),kstart,kend));
        colormap jet; 
%         set(gca,'FontSize',16); set(gcf,'position',[1077 230 412 378]);
%         shg; pause;
    end
    set(gca,'FontSize',16);
    set(gcf,'position',[1077 230 412 378]);
    title(['TE_{eff} = ' num2str(te_eff) 'ms, etl = ' num2str(etl)]);
    cbar=colorbar;
    cbar.Title.String = 'Echo #';
    xlabel('kx'); ylabel('ky');
    eval(['k_fill_image_' num2str(etl) '=k_fill_image;']);
        
    % Generate k-space
    brain_mese = brain_mese_etls{find(etls==etl),1};
    brain_kspace = zeros(dims);
    for ee=1:etl
%         kstart = (ee-1)*klines_per_echo+1;
%         kend   = kstart+klines_per_echo-1;
        kstart = k_echo_lines(klineorder(ee),1);
        kend = k_echo_lines(klineorder(ee),2);
        freqim = fftshift(fft2(brain_mese(:,:,ee)));
        brain_kspace(kstart:kend,:) = freqim(kstart:kend,:);
    end
    
    % Save simulated brain image
    brain_sim = abs(ifft2(fftshift(brain_kspace)));
    brain_sim_te80(:,:,tt) = brain_sim;
end


% Generate single-echo spin echo image w TE=80ms
brain_sese = zeros(dims);
for xx=1:dims(1)
    if mod(xx,10)==0, fprintf(' x=%d\n',xx); end
    for yy=1:dims(2)
        if mask(xx,yy)==1
            t1val = T1map(xx,yy);
            t2val = T2map(xx,yy);
            m0val = M0map(xx,yy);
            [Msig,Mss] = fsesignal_multi(t1val,t2val,80,TR,dfreq,1,m0val,1);
            brain_sese(xx,yy) = (Msig);
        end
    end
end

% Plot the images from (i) Single SE TE=80ms, (ii-v) FSE TE_eff=80ms with
%                                                       ETL={16,32,64,128}
f = figure('color','w','position',[158 503 1343 350]);
subplot(151);
imshow(brain_sese,[]);
title(['TE_{eff}=80ms, single-echo SE']);
set(gca,'FontSize',16);
subplot(152);
imshow(brain_sim_te80(:,:,1),[0.25 0.7]);
title(['TE_{eff}=80ms, etl=' num2str(etls(1))]);
set(gca,'FontSize',16);
subplot(153);
imshow(brain_sim_te80(:,:,2),[]);
title(['TE_{eff}=80ms, etl=' num2str(etls(2))]);
set(gca,'FontSize',16);
subplot(154);
imshow(brain_sim_te80(:,:,3),[]);
title(['TE_{eff}=80ms, etl=' num2str(etls(3))]);
set(gca,'FontSize',16);
subplot(155);
imshow(brain_sim_te80(:,:,4),[]);
title(['TE_{eff}=80ms, etl=' num2str(etls(4))]);
set(gca,'FontSize',16)
print(f,[plotdir filesep 'p2b-iv_brain-fse-sim-images.png'],'-dpng');

% Plot the kspace filling order from 
% (i) Single SE TE=80ms, (ii-v) FSE TE_eff=80ms with ETL={16,32,64,128}
f = figure('color','w','position',[8 617 1500 248]);
subplot(151);
imshow(ones(size(brain_sese)),[]);
title(['Single-echo SE']);
set(gca,'FontSize',12); set(gca,'colormap',jet);
cbar=colorbar; cbar.Title.String='Echo #';
ylabel('ky');xlabel('kx');
subplot(152);
imshow(k_fill_image_16,[1 16]);
title(['etl=' num2str(etls(1))]);
set(gca,'FontSize',12); set(gca,'colormap',jet);
cbar=colorbar; cbar.Title.String='Echo #';
ylabel('ky');xlabel('kx');
subplot(153);
imshow(k_fill_image_32,[1 32]);
title(['etl=' num2str(etls(2))]);
set(gca,'FontSize',12); set(gca,'colormap',jet);
cbar=colorbar; cbar.Title.String='Echo #';
ylabel('ky');xlabel('kx');
subplot(154);
imshow(k_fill_image_64,[1 64]);
title(['etl=' num2str(etls(3))]);
set(gca,'FontSize',12); set(gca,'colormap',jet);
cbar=colorbar; cbar.Title.String='Echo #';
ylabel('ky');xlabel('kx');
subplot(155);
imshow(k_fill_image_128,[1 128]);
title(['etl=' num2str(etls(4))]);
set(gca,'FontSize',12); set(gca,'colormap',jet);
cbar=colorbar; cbar.Title.String='Echo #';
ylabel('ky');xlabel('kx');

print(f,[plotdir filesep 'p2b-iv_brain-fse-sim-kspace-fill-order.png'],'-dpng');


%% Part c

% Set parameters to use
TR = 3000;
esp = 2;
numTRs = 1;
flipangle = pi;
etl = 128;
te_eff = 80;

brain_mese_etl128_esp2 = zeros([dims,etl]);
for xx=1:dims(1)
    if mod(xx,10)==0, fprintf(' x=%d\n',xx); end
    for yy=1:dims(2)
        if mask(xx,yy)==1
            t1val = T1map(xx,yy);
            t2val = T2map(xx,yy);
            m0val = M0map(xx,yy);
            [Msig,Mss] = fsesignal_multi(t1val,t2val,esp,TR,dfreq,etl,m0val,numTRs);
            brain_mese_etl128_esp2(xx,yy,:) = (Msig);
        end
    end
end


echotimes = (1:etl)*esp;
central_echo_ind = find(echotimes==te_eff);
klines_per_echo = dims(1)/etl;

klineorder = [104:128,1:39,40:103];

% Visualize k-space filling order:
k_fill_image = zeros(dims);
k_echo_lines = zeros(etl,2);
figure('position',[1005 251 484 357]); 
hold on; 
for ee=1:etl
    kstart = (ee-1)*klines_per_echo+1;
    kend   = kstart+klines_per_echo-1;
    k_fill_image(kstart:kend,:) = klineorder(ee);
    k_echo_lines(klineorder(ee),:) = [kstart,kend];
    imshow(k_fill_image,[0 etl]);
    title(sprintf('echo %d, ky=%d,%d',klineorder(ee),kstart,kend));
    colormap jet; 
%         set(gca,'FontSize',16); set(gcf,'position',[1077 230 412 378]);
%         shg; pause;
end
set(gca,'FontSize',16);
set(gcf,'position',[1077 230 412 378]);
title(['TE_{eff} = ' num2str(te_eff) 'ms, etl = ' num2str(etl)]);
cbar=colorbar;
cbar.Title.String = 'Echo #';
xlabel('kx'); ylabel('ky');
eval(['k_fill_image_' num2str(etl) '_esp2=k_fill_image;']);
    
% Generate k-space
brain_kspace = zeros(dims);
for ee=1:etl
%         kstart = (ee-1)*klines_per_echo+1;
%         kend   = kstart+klines_per_echo-1;
    kstart = k_echo_lines(klineorder(ee),1);
    kend = k_echo_lines(klineorder(ee),2);
    freqim = fftshift(fft2(brain_mese_etl128_esp2(:,:,ee)));
    brain_kspace(kstart:kend,:) = freqim(kstart:kend,:);
end

% Save simulated brain image
brain_sim = abs(ifft2(fftshift(brain_kspace)));
brain_sim_te80_etl128_esp2 = brain_sim;

% Plot comparing the (i) SE image (ii) FSE ETL=128 ESP=5ms image and 
% (iii) FSE ETL=128 ESP=2ms image 
% [notice the improved contrast in (iii) vs (ii)]
f = figure('color','w','position',[158 503 1343 350]);
subplot(131);
imshow(brain_sese,[]);
title(['TE_{eff}=80ms, single-echo SE']);
set(gca,'FontSize',16);
subplot(132);
imshow(brain_sim_te80(:,:,4),[]);
title(['TE_{eff}=80ms, etl=128, esp=5ms']);
set(gca,'FontSize',16);
subplot(133);
imshow(brain_sim_te80_etl128_esp2,[]);
title(['TE_{eff}=80ms, etl=128, esp=2ms']);
set(gca,'FontSize',16)
print(f,[plotdir filesep 'p2c_brain-fse-sim-images.png'],'-dpng');


% (i) Single SE TE=80ms, (ii-v) FSE TE_eff=80ms with ETL={16,32,64,128}
f = figure('color','w','position',[8 617 1500 248]);
subplot(131);
imshow(ones(size(brain_sese)),[]);
title(['Single-echo SE']);
set(gca,'FontSize',12); set(gca,'colormap',jet);
cbar=colorbar; cbar.Title.String='Echo #';
ylabel('ky');xlabel('kx');
subplot(132);
imshow(k_fill_image_128,[1 128]);
title(['TE_{eff}=80ms, etl=128, esp=5ms']);
set(gca,'FontSize',12); set(gca,'colormap',jet);
cbar=colorbar; cbar.Title.String='Echo #';
ylabel('ky');xlabel('kx');
subplot(133);
imshow(k_fill_image_128_esp2,[1 128]);
title(['TE_{eff}=80ms, etl=128, esp=2ms']);
set(gca,'FontSize',12); set(gca,'colormap',jet);
cbar=colorbar; cbar.Title.String='Echo #';
ylabel('ky');xlabel('kx');

print(f,[plotdir filesep 'p2c_brain-fse-sim-kspace-fill-order.png'],'-dpng');


