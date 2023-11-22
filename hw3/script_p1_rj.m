%% BME599 HW3 Problem 1: Partial Fourier Imaging

% Robert Jones
% 11-20-2023

clear
close all
% clc

writefigs = true;
figdir = '/Users/robertjones/Desktop/F23/599/hw/hw3/figures_p1_pocs_test';
if ~exist(figdir,'dir'), mkdir(figdir); end

addpath(genpath('/Users/robertjones/Desktop/F23/599/hw/hw3'));

%% Part (a): zero-filled reconstruction

%% load data + recon FS image
disp(' -Loading raw, fully-sampled, single-coil, complex kspace data...');
load('Data_Assignment3_Problem1.mat','kspaceData_SingleCoil');
ksp_fs = kspaceData_SingleCoil;
dims = size(ksp_fs);

disp(' -Fully-sampled FT image reconstruction...')
image_fs = ifftshift(ifft2(ifftshift(ksp_fs)));

f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(abs(image_fs),[0 3e-3]); title('FS (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(angle(image_fs),[]); title('FS (Phase)');
set(gca,'FontSize',20,'Colormap',hsv); colorbar;

fpng = [figdir filesep 'p1a_FS_image.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%% Undersampling for PF 5/8
disp(' -Undersampling for PF 5/8...');

DCoffset = ceil(dims(1)/2)+1;

PFrate = 5/8;
fprintf(' PF factor = %g\n',PFrate);

%% Make a PF sampling mask removing 3/8 of the kspace PE lines
[pelines,felines] = size(ksp_fs);
PF_start_line = 1;
PF_end_line = PFrate*pelines;
PFmask = zeros(dims);
PFmask(PF_start_line:PF_end_line,:) = 1;
PFlinesRecycle = PF_start_line:PF_end_line;

%% Make a POCS low-frequency mask (for phase estimation)
PocsMask = zeros(dims);
dcradius = (PFrate-0.5);
dclen = dcradius*dims(1);
dcstart = dims(2)*(0.5-dcradius);
dcend = dims(2)*(0.5+dcradius);
dcrange = dcstart+1:dcend;
filtlen = length(dcrange);
PocsMask(dcrange,:)=1;

% Plot the sampling masks
PlotMask = PocsMask+PFmask;
figure; 
imshow(PlotMask,[0 2]); colormap jet; colorbar;
title('PF mask==1, POCS mask==2, unsampled==0');

%% Generate PF kspace data (remove unsampled PE lines)

ksp_pf = ksp_fs .* PFmask;
fprintf(' nnz(PF kspace)=%g, numel(PF kspace)=%g, eff PF = %g\n', ...
    nnz(ksp_pf), numel(ksp_pf), ...
    nnz(ksp_pf)/numel(ksp_pf));

%% Zero-filled PF recon

disp(' -Zero-filled PF recon...');
% pf_zf_imageData_SingleCoil_Sym = fftshift(ifft2(ifftshift(pf_kspaceData_SingleCoil),'symmetric'));
image_pf_zp = ifftshift(ifft2(ifftshift(ksp_pf)));

f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(abs(image_pf_zp),[0 3e-3]); title('PF+ZP (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(angle(image_pf_zp),[]); title('PF+ZP (Phase)');
set(gca,'FontSize',20,'Colormap',hsv); colorbar;

fpng = [figdir filesep 'p1a_PF+ZF_image.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%% Compute difference between PF zero-filled and FullySampled image

diff_pf_zp = image_fs - image_pf_zp;

diff_pf_zp_magn = abs(diff_pf_zp);
diff_pf_zp_ph = angle(diff_pf_zp);

% diff_pf_zp_magn = abs(image_fs)-abs(image_pf_zp);
% diff_pf_zp_ph = angle(image_fs)-angle(image_pf_zp);


f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(diff_pf_zp_magn,[0 5e-4]); title('Diff. b/w FS & PF(5/8)+ZF (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow((diff_pf_zp_ph),[]); title('Diff. b/w FS & PF(5/8)+ZF (Phase)');
set(gca,'FontSize',20,'Colormap',gray); colorbar;

fpng = [figdir filesep 'p1a_PF+ZF_difference.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

%% Part (b): POCS Conjugate Phase Reconstruction

[ pocs_recon_image, pocs_recon_kspace ] = pocs_recon( ksp_pf, 5/8, 25, 0 );

f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(abs(pocs_recon_image),[0 3.5e-3]); title('PF+POCS (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
% imshow(imag(pocs_recon_image),[]); title('Pocs recon (Phase)');
imshow(angle(pocs_recon_image),[]); title('PF+POCS (Phase)');
set(gca,'FontSize',20,'Colormap',hsv); colorbar;

fpng = [figdir filesep 'p1a_PF+POCS-recon_image.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%% Compute difference between PF zero-filled and FullySampled image
diff_pf_pocs = image_fs - pocs_recon_image;
diff_pf_pocs_magn = abs(diff_pf_pocs);
diff_pf_pocs_ph = angle(diff_pf_pocs);

% diff_pf_pocs_magn = abs(image_fs)-abs(pocs_recon_image);
% diff_pf_pocs_ph = angle(image_fs)-angle(pocs_recon_image);


f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(diff_pf_pocs_magn,[0 5e-4]); title('Diff. b/w FS & PF(5/8)+POCS (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(diff_pf_pocs_ph,[]); title('Diff. b/w FS & PF(5/8)+POCS (Phase)');
set(gca,'FontSize',20,'Colormap',gray); colorbar;
% print(f,'p1b-difference_pocs-vs-fs.png'],'-dpng','-r300');

fpng = [figdir filesep 'p1a_PF+POCS_difference.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%% Difference bw ZP and FS
diff_pf_zp = image_fs - image_pf_zp;

diff_pf_zp_magn = abs(diff_pf_zp);
diff_pf_zp_ph = angle(diff_pf_zp);

% diff_pf_zp_magn = abs(image_fs)-abs(image_pf_zp);
% diff_pf_zp_ph = angle(image_fs)-angle(image_pf_zp);


f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(diff_pf_zp_magn,[0 5e-4]); title('Diff. b/w FS & PF(5/8)+zP (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(diff_pf_zp_ph,[]); title('Diff. b/w FS & PF(5/8)+ZP (Phase)');
set(gca,'FontSize',20,'Colormap',gray); colorbar;
% print(f,'p1b-difference_pocs-vs-fs.png'],'-dpng','-r300');

fpng = [figdir filesep 'p1a_PF+POCS_difference.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%%% Compute RMSEs of ZP vs POCS

rmse.magn.zp = sqrt(mean( (abs(image_fs(:))-abs(image_pf_zp(:))).^2 ) );
rmse.magn.pocs = sqrt(mean( (abs(image_fs(:))-abs(pocs_recon_image(:))).^2 ) );

rmse.ph.zp = sqrt(mean( (angle(image_fs(:))-angle(image_pf_zp(:))).^2 ) );
rmse.ph.pocs = sqrt(mean( (angle(image_fs(:))-angle(pocs_recon_image(:))).^2 ) );




% %% OLD REMOVE:!!!
% 
% % Make hanning filter for PE dir with width == reference DC region
% hfilt = hann(filtlen);
% hfiltpad = zeros(dims(1),1);
% hfiltpad(dcrange) = hfilt;
% HanningWindow = repmat(hfiltpad,1,dims(2));
% % figure; imshow(HanningWindow);
% 
% %%%% reduce width of hanning filter to half reference DC region
% % hfilt = hann(25);
% % hfiltpad = zeros(dims(1),1);
% % hfiltpad(88:112) = hfilt;
% % HanningWindow = repmat(hfiltpad,1,dims(2));
% % figure; imshow(HanningWindow);
% 
% % now make hanning window for FE (fully sampled) dim:
% hfilt_freq = reshape(hann(dims(2)),1,[]);
% HanningWindow_freq = repmat(hfilt_freq,dims(1),1);
% % figure; imshow(HanningWindow_freq);
% 
% HanningWindow2D = HanningWindow .* HanningWindow_freq;
% figure; imshow(HanningWindow2D); title('2d hanning window');
% 
% 
% ksp_pf_filt = ksp_pf .* HanningWindow2D; %HanningWindow;
% image_pf_filt = ifftshift(ifft2(ifftshift(ksp_pf_filt)));
% 
% 
% 
% 
% 
% figure; imshow(abs(ksp_pf_filt),[0 0.0001]); colormap jet;
% 
% % PocsEstimatedPhase = imag(imagePartialFilt);
% PocsEstimatedPhase = angle(image_pf_filt);
% 
% f = figure('position', [163 86 429 780],'color','w');
% subplot(211);
% imshow(abs(image_pf_filt),[]); title('POCS Ph Ref Hann (Magn)');
% set(gca,'FontSize',20); colorbar;
% subplot(212);
% imshow(angle(image_pf_filt),[]); title('POCS Ph Ref Hann (Phase)');
% set(gca,'FontSize',20,'Colormap',hsv); colorbar;
% 
% fpng = [figdir filesep 'p1a_POCS-ref_image.png'];
% if ~exist(fpng,'file') || writefigs
%     print(f,fpng,'-dpng');
% end
% 
% niters = 25;
% itererrs = zeros(niters,2);
% gt1 = reshape(abs(image_fs),[],1);
% % gt2 = reshape(imag(imageData_SingleCoil),[],1);
% gt2 = reshape(angle(image_fs),[],1);
% itermaps = zeros([dims 2 niters]);
% 
% tmpkspace = ksp_pf;
% % keeplines = 1:PFrate*dims(1); replacelines = PFrate*dims(1)+1:dims(1);
% for nn=1:niters
%     fprintf(' POCS iter %d\n',nn);
%     tmpimage = ifftshift(ifft2(ifftshift(tmpkspace)));
%     itermaps(:,:,1,nn) = abs(tmpimage);
%     itermaps(:,:,2,nn) = imag(tmpimage);
% 
%     e1 = reshape(abs(tmpimage),[],1);
%     e2 = reshape(angle(tmpimage),[],1);
% 
%     itererrs(nn,1) = sqrt(mean((gt1-e1).^2));
%     itererrs(nn,2) = sqrt(mean((gt2-e2).^2));
%     
% %     tmpimage = real(tmpimage) + 1i*PocsEstimatedPhase;
%     tmpimage = abs(tmpimage) .* exp(1i.*PocsEstimatedPhase);
%     e1 = reshape(abs(tmpimage),[],1);
%     e2 = reshape(exp(1i.*angle(tmpimage)),[],1);
% 
%     itererrs(nn,1) = sqrt(mean((gt1-e1).^2));
%     itererrs(nn,2) = sqrt(mean((gt2-e2).^2));
% 
%     tmpkspace = fftshift(fft2(fftshift(tmpimage)));
%     tmpkspace(PFlinesRecycle,:) = ksp_fs(PFlinesRecycle,:);
% end
% 
% figure; 
% subplot(121); plot(itererrs(:,1))
% subplot(122); plot(real(itererrs(:,2)))
% 
% 
% %% POCS image reconstruction
% 
% pocs_recon_image = ifftshift(ifft2(ifftshift(tmpkspace)));
% 
% f = figure('position', [163 86 429 780],'color','w');
% subplot(211);
% imshow(abs(pocs_recon_image),[0 3.5e-3]); title('PF+POCS (Magn)');
% set(gca,'FontSize',20); colorbar;
% subplot(212);
% % imshow(imag(pocs_recon_image),[]); title('Pocs recon (Phase)');
% imshow(angle(pocs_recon_image),[]); title('PF+POCS (Phase)');
% set(gca,'FontSize',20,'Colormap',hsv); colorbar;
% 
% fpng = [figdir filesep 'p1a_PF+POCS-recon_image.png'];
% if ~exist(fpng,'file') || writefigs
%     print(f,fpng,'-dpng');
% end
% 
% % 
% % [otherpocs_image, otherpocs_kspFull] = pocs( ksp_pf, 25, 1 );
% % 
% % f = figure('position', [163 86 429 780],'color','w');
% % subplot(211);
% % imshow(abs(otherpocs_image),[0 3.5e-3]); title('PF+POCS(other) (Magn)');
% % set(gca,'FontSize',20); colorbar;
% % subplot(212);
% % % imshow(imag(pocs_recon_image),[]); title('Pocs recon (Phase)');
% % imshow(angle(otherpocs_image),[]); title('PF+POCS(other) (Phase)');
% % set(gca,'FontSize',20,'Colormap',hsv); colorbar;
% % 
% % fpng = [figdir filesep 'p1a_PF+otherPOCS-recon_image.png'];
% % if ~exist(fpng,'file') || writefigs
% %     print(f,fpng,'-dpng');
% % end
% 
% 
% %% Compute difference between PF zero-filled and FullySampled image
% diff_pf_pocs = image_fs - pocs_recon_image;
% 
% diff_pf_pocs_magn = abs(diff_pf_pocs);
% diff_pf_pocs_ph = angle(diff_pf_pocs);
% 
% % diff_pf_pocs_magn = abs(image_fs)-abs(pocs_recon_image);
% % diff_pf_pocs_ph = angle(image_fs)-angle(pocs_recon_image);
% 
% 
% f = figure('position', [163 86 429 780],'color','w');
% subplot(211);
% imshow(diff_pf_pocs_magn,[0 5e-4]); title('Diff. b/w FS & PF(5/8)+POCS (Magn)');
% set(gca,'FontSize',20); colorbar;
% subplot(212);
% imshow(diff_pf_pocs_ph,[]); title('Diff. b/w FS & PF(5/8)+POCS (Phase)');
% set(gca,'FontSize',20,'Colormap',gray); colorbar;
% % print(f,'p1b-difference_pocs-vs-fs.png'],'-dpng','-r300');
% 
% fpng = [figdir filesep 'p1a_PF+POCS_difference.png'];
% if ~exist(fpng,'file') || writefigs
%     print(f,fpng,'-dpng');
% end
% 
% 
% %% Difference b/w otherpocs + fs
% diff_pf_pocsother = image_fs - otherpocs_image;
% 
% diff_pf_pocsother_magn = abs(diff_pf_pocsother);
% diff_pf_pocsother_ph = angle(diff_pf_pocsother);
% 
% % diff_pf_pocs_magn = abs(image_fs)-abs(pocs_recon_image);
% % diff_pf_pocs_ph = angle(image_fs)-angle(pocs_recon_image);
% 
% 
% f = figure('position', [163 86 429 780],'color','w');
% subplot(211);
% imshow(diff_pf_pocsother_magn,[0 5e-4]); title('Diff. b/w FS & PF(5/8)+POCS* (Magn)');
% set(gca,'FontSize',20); colorbar;
% subplot(212);
% imshow(diff_pf_pocsother_ph,[]); title('Diff. b/w FS & PF(5/8)+POCS* (Phase)');
% set(gca,'FontSize',20,'Colormap',gray); colorbar;
% % print(f,'p1b-difference_pocs-vs-fs.png'],'-dpng','-r300');
% 
% fpng = [figdir filesep 'p1a_PF+POCSother_difference.png'];
% if ~exist(fpng,'file') || writefigs
%     print(f,fpng,'-dpng');
% end
% 
% 
% %%% Compute RMSEs of ZP vs POCS
% 
% rmse.magn.zp = sqrt(mean( (abs(image_fs(:))-abs(image_pf_zp(:))).^2 ) );
% rmse.magn.pocs = sqrt(mean( (abs(image_fs(:))-abs(pocs_recon_image(:))).^2 ) );
% 
% rmse.ph.zp = sqrt(mean( (angle(image_fs(:))-angle(image_pf_zp(:))).^2 ) );
% rmse.ph.pocs = sqrt(mean( (angle(image_fs(:))-angle(pocs_recon_image(:))).^2 ) );
% 
% 
% % diff_pf_zp = image_fs - image_pf_zp;
% % diff_pf_zp_magn = abs(image_fs) - abs(image_pf_zp);
% % diff_pf_zp_ph = angle(image_fs) - angle(image_pf_zp);
% % 
% % f = figure('position', [245 490 1115 376],'color','w');
% % subplot(121);
% % imshow(diff_pf_zp_magn,[]); title('Diff. b/w FS & PF(5/8)+ZF (Magn)');
% % set(gca,'FontSize',20); colorbar;
% % subplot(122);
% % imshow((diff_pf_zp_ph),[]); title('Diff. b/w FS & PF(5/8)+ZF (Phase)');
% % set(gca,'FontSize',20,'Colormap',jet); colorbar;
% % print(f,'p1b-difference_pocs-vs-fs.png'],'-dpng','-r300');
% 

