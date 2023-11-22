%% BME599 HW3 Problem 2: SENSE

% Robert Jones
% 11-20-2023

clear
close all
% clc

writefigs = true;
figdir = '/Users/robertjones/Desktop/F23/599/hw/hw3/figures_p2_sense_FINAL2';
if ~exist(figdir,'dir'), mkdir(figdir); end

addpath(genpath('/Users/robertjones/Desktop/F23/599/hw/hw3'));

%% Part (a): Fully-sampled image

%% Load data
disp(' -Loading raw, fully-sampled data...');
load('Data_Assignment3_Problem2.mat');
ksp_fs = kspaceData;
dims = size(ksp_fs,1:2);
ncoils = size(ksp_fs,3);

% Reconstruct coil images
disp(' -Fully-sampled FT image reconstruction...')
coilimages_fs = zeros(size(ksp_fs)); %ifftshift(ifft2(ifftshift(ksp_fs)));
for nn=1:ncoils
    coilimages_fs(:,:,nn) = ifftshift(ifft2(ifftshift(ksp_fs(:,:,nn))));
end

% Plot the coil sens maps and FS coil images
coiltile = mat2gray(abs(coilimages_fs),[0 3e-3]);
it1 = imtile(coiltile,'GridSize', [2 4]);
senstile = mat2gray(abs(coilmaps),[0 1]);
it2 = imtile(senstile,'GridSize', [2 4]);

f = figure('color','w', 'position',[476 206 587 660]); %[15 270 650 596]);
subplot(211); 
imshow(it2); title('coil sens maps'); set(gca,'FontSize',24);
subplot(212); 
imshow(it1); title('coil images'); set(gca,'FontSize',24);
fpng = [figdir filesep 'p2a_tiled_coilmaps+coilimages.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

%% Coil combination FS data

% Do coil combination using coil images + coil sens maps
% image_combine_fs_ssa = sqrt(sum(abs(coilimages_fs).^2,3));
image_combine_fs = (sum(conj(coilmaps).*coilimages_fs,3) ./ sum(conj(coilmaps).*coilmaps,3));

% Plot the coil combined image magn + ph
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow((abs(image_combine_fs)),[0 4e-3]); title('FS coil combined (Magn)');
% imshow((abs(image_combine_fs)),[]); title('FS coil combined (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(angle(image_combine_fs),[]); title('FS coil combined (Phase)');
set(gca,'FontSize',20,'Colormap',hsv); colorbar;

fpng = [figdir filesep 'p2a_FS-coil-combined-image.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%% Undersampling R=2 ( along PE dir )

Rfactor = 2;

% undersample the kpsace data by R=2
ksp_r2 = ksp_fs;
ksp_r2(2:Rfactor:end,:,:)=0;
max(abs(ksp_r2(101,:)))

% Recon undersampled R=2 coil images
coilimages_r2 = zeros(size(ksp_r2)); 
for nn=1:ncoils
    coilimages_r2(:,:,nn) = ifftshift(ifft2(ifftshift(ksp_r2(:,:,nn))));
end

% Combine undersampled R=2 coil images
image_combine_r2 = sum(conj(coilmaps).*coilimages_r2,3) ./ sum(conj(coilmaps).*coilmaps,3);
% image_combine_r2 = sqrt(sum(abs(coilimages_r2).^2,3));

% Plot the R=2 undersampled (direct recon) imgaes
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(abs(image_combine_r2),[]); title('R=2 coil combined (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(angle(image_combine_r2),[]); title('R=2 coil combined (Phase)');
set(gca,'FontSize',20,'Colormap',hsv); colorbar;

fpng = [figdir filesep 'p2b_R=2_coil-combined-image.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%% Implement SENSE R=2

alias_px = dims(1)/Rfactor; % = round(dims(1)/Rfactor;
alias_offset = alias_px/2;

% Crop the undersampled coil images to the "Nyquist" FOV
coilimages_r2_crop = (coilimages_r2(alias_offset+1:dims(1)-alias_offset,:,:));
cropdims = size(coilimages_r2_crop);

sense_recon_image.r2 = zeros(dims);

for xx=1:cropdims(1)
    fprintf(' x=%d\n',xx);
    linenum = alias_offset+xx;
    x1 = linenum;
    x2 = x1+alias_px;
    if x2>dims(1), x2=x2-dims(1); end
    for yy=1:cropdims(2)
        clear sensecoef imagecoef
        imagecoef = squeeze(coilimages_r2_crop(xx,yy,:));
        sensecoef = zeros(ncoils,Rfactor);
        sensecoef(:,1) = coilmaps(x1,yy,:);
        sensecoef(:,2) = coilmaps(x2,yy,:);       
        rho = pinv(sensecoef)*imagecoef;
        % (Can also use: lsqminnorm(sensecoef,imagecoef);)
        sense_recon_image.r2(x1,yy) = rho(1);
        sense_recon_image.r2(x2,yy) = rho(2);
    end   
end

%%%% You can use this to manually unwrap the estimated rhos to full FOV :
%    (This is equivalent to sense_recon_image)
% sense_recon_image = cat(1, abs(rhos(51:100,:,2)), abs(rhos(:,:,1)), abs(rhos(1:50,:,2)) );

% Scale by 2 to account for R=2 unaliasing
sense_recon_image.r2 = sense_recon_image.r2*Rfactor;
sense_recon_image_magn.r2 = abs(sense_recon_image.r2);
sense_recon_image_ph.r2 = angle(sense_recon_image.r2);

% Plot the sense recon image magn and phase
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(sense_recon_image_magn.r2,[]); 
title(sprintf('R=%d SENSE recon magn',Rfactor));colorbar;
set(gca,'FontSize',20);
subplot(212);
imshow(abs(image_combine_r2),[]); 
title(sprintf('R=%d Direct recon magn',Rfactor));
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2_R=' num2str(Rfactor) '_magn-images_SENSE+Direct.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

% Plot the sense recon image magn vs FS image magn
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(sense_recon_image_magn.r2,[]); 
title(sprintf('R=%d SENSE recon magn',Rfactor));colorbar;
set(gca,'FontSize',20);
subplot(212);
imshow(abs(image_combine_fs),[]); 
title(sprintf('FS recon magn'));
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2_R=' num2str(Rfactor) '_magn-images_SENSE+FS.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(sense_recon_image_magn.r2,[]); 
title(sprintf('R=%d SENSE recon magn',Rfactor));colorbar;
set(gca,'FontSize',20);
subplot(212);
imshow(abs(image_combine_r2),[]); 
title(sprintf('R=%d Direct recon magn',Rfactor));
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2_R=' num2str(Rfactor) '_magn-images_SENSE+Direct.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end




%% Compute magnitude difference between SENSE and FS recon images
% complex difference
diff_r2.sense = image_combine_fs - sense_recon_image.r2;
% magnitude of complex difference
diff_r2.sense_magn = abs(diff_r2.sense);
% phase angle of complex difference
diff_r2.sense_ph = angle(diff_r2.sense);

%% Compute magnitude difference between direct(ZF) and FS recon images
% complex difference
diff_r2.direct = image_combine_fs - image_combine_r2;
% magnitude of complex difference
diff_r2.direct_magn = abs(diff_r2.direct);
% phase angle of complex difference
diff_r2.direct_ph = angle(diff_r2.direct);

%% Plot magnitude differences

f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(diff_r2.sense_magn,[]); 
title(sprintf('Magn. Diff., R=%d, SENSE',Rfactor));
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(diff_r2.direct_magn,[]); 
title(sprintf('Magn. Diff., R=%d, Direct',Rfactor));
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2_R=' num2str(Rfactor) '_magn-diff_SENSE+Direct-vs-FS.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%% SENSE R=4

% Set acc factor
Rfactor = 4;

% Set reference region cropping
r4_pe_crop = 76:125;
pe_start = 1;

ksp_r4 = zeros(size(ksp_fs));
ksp_r4(pe_start:Rfactor:end,:,:)=ksp_fs(pe_start:Rfactor:end,:,:);
% max(abs(ksp_r4(101,:)));

% Recon undersampled R=2 coil images
coilimages_r4 = zeros(size(ksp_r2)); 
for nn=1:ncoils
    coilimages_r4(:,:,nn) = ifftshift(ifft2(ifftshift(ksp_r4(:,:,nn))));
end
% Combine undersampled R=2 coil images
% temp = zeros(dims);
% for nn=1:ncoils
%     temp = temp + conj(coilmaps(:,:,nn)).*coilimages_r4(:,:,nn);
% end
% image_combine_r4 = temp.^(1/2);
image_combine_r4 = sum(conj(coilmaps).*coilimages_r4,3) ./ sum(conj(coilmaps).*coilmaps,3);

% Plot the r=4 coil combined image magn + phase
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(abs(image_combine_r4),[]); title('R=4 coil combined (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(angle(image_combine_r4),[]); title('R=4 coil combined (Phase)');
set(gca,'FontSize',20,'Colormap',hsv); colorbar;

fpng = [figdir filesep 'p2d_R=4_coil-combined-image.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

% Plot the r=4 coil combined image magn + cropped magn
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(abs(image_combine_r4),[]); title('R=4 coil combined (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(abs(image_combine_r4(r4_pe_crop,:)),[]); title('R=4 coil combined (Magn zoom)');
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2d_R=4_coil-combined-image-magn+zoom.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

% Plot tiles of coil images r=4
coiltile_r4 = mat2gray(abs(coilimages_r4),[0 3e-3]);
it1_r4 = imtile(coiltile_r4,'GridSize', [2 4]);
it1_r4_crop = imtile(coiltile_r4(r4_pe_crop,:,:),'GridSize', [2 4]);

f = figure('color','w', 'position',[476 206 587 660]); %[15 270 650 596]);
subplot(211); 
imshow(it2); title('coil sens maps'); set(gca,'FontSize',24);
subplot(212); 
imshow(it1_r4); title('coil images R=4'); set(gca,'FontSize',24);

fpng = [figdir filesep 'p2d_tiled_coilmaps+coilimages_R=4.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

f = figure('color','w', 'position',[476 206 587 660]); %[15 270 650 596]);
imshow(it1_r4_crop); title('coil images R=4'); set(gca,'FontSize',24);
set(f,'Position',[476 206 695 250]);
fpng = [figdir filesep 'p2d_tiled_coilimages_R=4_zoom.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end




%% Implement SENSE

% Crop the coil images to the Nyquist FOV
coilimages_r4_crop = (coilimages_r4(r4_pe_crop,:,:));
r4_fov = length(r4_pe_crop);
cropdims_r4 = size(coilimages_r4_crop);

% rhos = zeros(cropdims_r4(1),cropdims_r4(2),Rfactor);
sense_recon_image.r4 = zeros(dims);

% Iterative SENSE recon: loop thru PE lines, for each loop thru FE samples
for xx=1:cropdims_r4(1)
    % The current undersampled PE line
    fprintf(' x=%d\n',xx);
    % Get theoretical unaliased (full FOV) PE line numbers
    linenum = r4_pe_crop(1)-1+xx;
    xv = mod([linenum, linenum+r4_fov, linenum+2*r4_fov, linenum+3*r4_fov],dims(1));
    xv(xv==0) = dims(1);
    [xv,xvind] = sort(xv); % Sort to avoid indexing problems :)
    for yy=1:cropdims(2)
        % The current FE step
        clear sensecoef imagecoef
        % Get coil image intensities
        imagecoef = squeeze(coilimages_r4_crop(xx,yy,:));
        % Get coil sens map intensities
        sensecoef = zeros(ncoils,Rfactor);
        sensecoef(:,1) = coilmaps(xv(1),yy,:);
        sensecoef(:,2) = coilmaps(xv(2),yy,:);
        sensecoef(:,3) = coilmaps(xv(3),yy,:);
        sensecoef(:,4) = coilmaps(xv(4),yy,:);
%         sensecoef = reshape(squeeze(coilmaps(xv,yy,:)),ncoils,Rfactor);
        % Estimate unaliased image intensities
        rho = pinv(sensecoef)*imagecoef;
        % Store resulting unaliased signal intensity estimations
        sense_recon_image.r4(xv(1),yy) = rho(1);
        sense_recon_image.r4(xv(2),yy) = rho(2);
        sense_recon_image.r4(xv(3),yy) = rho(3);
        sense_recon_image.r4(xv(4),yy) = rho(4);
%         sense_recon_image.r4(xv,yy) = rho;
    end   
end

sense_recon_image.r4 = sense_recon_image.r4 * Rfactor;
sense_recon_image_magn.r4 = abs(sense_recon_image.r4);
sense_recon_image_ph.r4 = angle(sense_recon_image.r4);

% Plot the sense recon image magn and phase
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(sense_recon_image_magn.r4,[]); 
title(sprintf('R=%d SENSE recon magn',Rfactor));colorbar;
set(gca,'FontSize',20);
subplot(212);
imshow(abs(image_combine_r4),[]); 
title(sprintf('R=%d Direct recon magn',Rfactor));
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2d_R=4_magn-images_SENSE+Direct.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

% Plot the sense recon image magn vs FS image magn 
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(sense_recon_image_magn.r4,[]); 
title(sprintf('R=%d SENSE recon magn',Rfactor));colorbar;
set(gca,'FontSize',20);
subplot(212);
imshow(abs(image_combine_fs),[]); 
title(sprintf('FS recon magn'));
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2d_R=4_magn-images_SENSE+FS.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end



f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(sense_recon_image_magn.r4,[]); 
title(sprintf('R=%d SENSE recon magn',Rfactor));colorbar;
set(gca,'FontSize',20);
subplot(212);
imshow(abs(image_combine_r4),[]); 
title(sprintf('R=%d Direct recon magn',Rfactor));
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2d_R=4_magn-images_SENSE+Direct.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end


%% Compute magnitude difference between SENSE and FS recon images
% complex difference
diff_r4.sense = image_combine_fs - (sense_recon_image.r4);
% magnitude of complex difference
diff_r4.sense_magn = abs(diff_r4.sense);
% phase angle of complex difference
diff_r4.sense_ph = angle(diff_r4.sense);

%% Compute magnitude difference between direct(ZF) and FS recon images
% complex difference
diff_r4.direct = image_combine_fs - image_combine_r4;
% magnitude of complex difference
diff_r4.direct_magn = abs(diff_r4.direct);
% phase angle of complex difference
diff_r4.direct_ph = angle(diff_r4.direct);

%% Plot magnitude differences

f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(diff_r4.sense_magn,[]); 
title(sprintf('Magn. Diff., R=%d, SENSE',Rfactor));
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(diff_r4.direct_magn,[]); 
title(sprintf('Magn. Diff., R=%d, Direct',Rfactor));
set(gca,'FontSize',20); colorbar;

fpng = [figdir filesep 'p2d_R=4_magn-diff_SENSE+Direct-vs-FS.png'];
if ~exist(fpng,'file') || writefigs
    print(f,fpng,'-dpng');
end

