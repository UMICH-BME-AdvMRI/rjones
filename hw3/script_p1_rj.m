%% BME599 HW3 Problem 1: Partial Fourier Imaging

% Robert Jones
% 11-20-2023

clear
close all
% clc

writefigs = true;
figdir = '/Users/robertjones/Desktop/F23/599/hw/hw3/figures_p1_pocs_FINAL2';
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


%% Undersampling for PF 5/8
disp(' -Undersampling for PF 5/8...');

DCoffset = ceil(dims(1)/2)+1;
PFrate = 5/8;
fprintf(' PF factor = %g\n',PFrate);

%% Make a PF sampling mask removing 3/8 of the kspace PE lines
[ PFmask ] = makePartialFourierSamplingMask( dims, PFrate );

%% Make a POCS low-frequency mask (for phase estimation)
[ PocsMask ] = makePocsReferenceMask( dims, PFrate );

%% Generate PF kspace data (remove unsampled PE lines)
ksp_pf = ksp_fs .* PFmask;
fprintf(' nnz(PF kspace)=%g, numel(PF kspace)=%g, eff PF = %g\n', ...
    nnz(ksp_pf), numel(ksp_pf), ...
    nnz(ksp_pf)/numel(ksp_pf));

%% Zero-filled PF recon
disp(' -Zero-filled PF recon...');
image_pf_zp = ifftshift(ifft2(ifftshift(ksp_pf)));



%% Compute difference between PF zero-filled and FullySampled image
    
if 1 == 1
    % the difference image
    diff_pf_zp = image_fs - image_pf_zp;
    % the magn and phase of the difference image
    diff_pf_zp_magn = abs(diff_pf_zp);
    diff_pf_zp_ph = angle(diff_pf_zp);
else
    % alternatively, take magn and phase and then compute difference
    diff_pf_zp_magn = abs(image_fs)-abs(image_pf_zp);
    diff_pf_zp_ph = angle(image_fs)-angle(image_pf_zp);
end


f = figure('position', [163 165 1299 701],'color','w');

subplot(2,3,1);
imshow(abs(image_fs),[0 3e-3]); title({'FS Magn.'});
set(gca,'FontSize',18); colorbar;
subplot(2,3,4);
imshow(angle(image_fs),[]); title({'FS Phase'});
set(gca,'FontSize',18,'Colormap',hsv); colorbar;

subplot(2,3,2);
imshow(abs(image_pf_zp),[0 3e-3]); title({'PF[5/8]+ZF Magn.'});
set(gca,'FontSize',18); colorbar;
subplot(2,3,5);
imshow(angle(image_pf_zp),[]); title({'PF[5/8]+ZF Phase'});
set(gca,'FontSize',18,'Colormap',hsv); colorbar;

subplot(2,3,3);
imshow(diff_pf_zp_magn,[0 5e-4]); title({'Magn. Difference','FS vs PF[5/8]+ZF'});
set(gca,'FontSize',18); colorbar;
subplot(2,3,6);
imshow((diff_pf_zp_ph),[]); title({'Phase Difference','FS vs PF[5/8]+ZF'});
set(gca,'FontSize',18,'Colormap',hsv); colorbar;

fpng = [figdir filesep 'p1a_images.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end



%% Part (b): POCS Conjugate Phase Reconstruction
niters_pocs = 500;
fprintf(' Running POCS reconstruction with %d iters...\n',niters_pocs);
[ pocs_recon_image, pocs_recon_kspace ] = pocs_recon( ksp_pf, PFrate, niters_pocs, 0 );

%% Compute difference between PF zero-filled and FullySampled image
if 1 == 1
    % the difference image
    diff_pf_pocs = image_fs - pocs_recon_image;
    % the magn and phase of the difference image
    diff_pf_pocs_magn = abs(diff_pf_pocs);
    diff_pf_pocs_ph = angle(diff_pf_pocs);
else
    % alternatively, take magn and phase and then compute difference
    diff_pf_pocs_magn = abs(image_fs)-abs(diff_pf_pocs);
    diff_pf_pocs_ph = angle(image_fs)-angle(diff_pf_pocs);
end


f = figure('position', [163 165 1299 701],'color','w');

subplot(2,3,1);
imshow(abs(image_fs),[0 3e-3]); title({'FS Magn.'});
set(gca,'FontSize',18); colorbar;
subplot(2,3,4);
imshow(angle(image_fs),[]); title({'FS Phase'});
set(gca,'FontSize',18,'Colormap',hsv); colorbar;

subplot(2,3,2);
imshow(abs(pocs_recon_image),[0 3e-3]); title({'PF[5/8]+POCS Magn.'});
set(gca,'FontSize',18); colorbar;
subplot(2,3,5);
imshow(angle(pocs_recon_image),[]); title({'PF[5/8]+POCS Phase'});
set(gca,'FontSize',18,'Colormap',hsv); colorbar;

subplot(2,3,3);
imshow(diff_pf_pocs_magn,[0 5e-4]); title({'Magn. Difference','FS vs PF[5/8]+POCS'});
set(gca,'FontSize',18); colorbar;
subplot(2,3,6);
imshow((diff_pf_pocs_ph),[]); title({'Phase Difference','FS vs PF[5/8]+POCS'});
set(gca,'FontSize',18,'Colormap',hsv); colorbar;

fpng = [figdir filesep 'p1b_images_' num2str(niters_pocs) 'iters.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

%%% Compute RMSEs of ZP vs POCS
brainmask = abs(image_fs)>4e-4;

rmse.magn.zp = sqrt(mean( (abs(image_fs(brainmask))-abs(image_pf_zp(brainmask))).^2 ) );
rmse.magn.pocs = sqrt(mean( (abs(image_fs(brainmask))-abs(pocs_recon_image(brainmask))).^2 ) );
nrmse.magn.zp = sqrt(mean( ((abs(image_fs(brainmask))-abs(image_pf_zp(brainmask)))./(abs(image_fs(brainmask)))).^2 ) );
nrmse.magn.pocs = sqrt(mean( ((abs(image_fs(brainmask))-abs(pocs_recon_image(brainmask)))./(abs(image_fs(brainmask)))).^2 ) );
rmse.ph.zp = sqrt(mean( (angle(image_fs(brainmask))-angle(image_pf_zp(brainmask))).^2 ) );
rmse.ph.pocs = sqrt(mean( (angle(image_fs(brainmask))-angle(pocs_recon_image(brainmask))).^2 ) );

