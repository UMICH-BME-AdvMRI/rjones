function [ pocs_image, pocs_kspace ] = pocs_recon( kspace, pf, niters, makeplots )
% [ pocs_image, pocs_kspace ] = pocs_recon( kspace, pf, niters, makeplots )
% 
% set showprogress=false to not make plots
% 

showprogress = true; pause off; plotpause = 0.01;
useconj = true;

if nargin<1
    fprintf(' USAGE:\n');
end
dims = size(kspace);

if nargin<4
    makeplots = false;
else
    makeplots = logical(makeplots);
end
if nargin<3
    niters = 25;
end
if nargin<2
    pf = nnz(find(max(abs(kspace),[],2)>0))/dims(1);
end

kprec = class(kspace);

dc_offset = floor(dims(1)/2)+1;

% Find fully-sampled region for phase estimation
nreflines = 2*round( (pf-0.5)*dims(1) );
reflines = dc_offset-(nreflines/2) : dc_offset+(nreflines/2)-1;

% Find acquired PE lines
sampledPeLines = 1:reflines(end);

% Make hanning filter for PE dir with width == reference DC region
hfilt_pe = zeros(dims(1),1);
hfilt_pe(reflines) = hann(nreflines);
hfilt_pe = repmat(hfilt_pe,1,dims(2));

% Now make hanning window for FE (fully sampled) dim:
hfilt_fe = reshape(hann(dims(2)),1,[]);
hfilt_fe = repmat(hfilt_fe,dims(1),1);

% Combine to get 2D hanning filter
hannfilt = hfilt_pe .* hfilt_fe;
if makeplots,figure;imshow(hannfilt);title('2D hanning filt');drawnow;end

% Generate (zp) symmetric reference kspace
kspace_filt = kspace .* hannfilt; %HanningWindow;
% Recon reference image
image_pf_filt = ifftshift(ifft2(ifftshift(kspace_filt)));

if makeplots,figure;imshow(abs(kspace_filt)>0);colormap jet;end

if makeplots
f = figure('position', [163 86 429 780],'color','w');
subplot(211);
imshow(abs(image_pf_filt),[]); title('POCS Ph Ref Hann (Magn)');
set(gca,'FontSize',20); colorbar;
subplot(212);
imshow(angle(image_pf_filt),[]); title('POCS Ph Ref Hann (Phase)');
set(gca,'FontSize',20,'Colormap',hsv); colorbar;
% 
% fpng = 'p1a_POCS-ref_image.png';
% if ~exist(fpng,'file') || writefigs
%     print(f,fpng,'-dpng');
% end
end

% extract reference (low-res) phase for pocs iterations
ph = exp(1i*angle(image_pf_filt));
if 1 == 1 %useconj
ph = ph./prod(dims);
end

% Initialize image with magn from ZP PF FT recon and phase from reference
if useconj
iminit = fftshift(fft2(fftshift( conj(kspace) )));
else
iminit = fftshift(fft2(fftshift( (kspace) )));
end
iminit = abs(iminit) .* ph;

tmpim = iminit;

if showprogress
    f=figure; drawnow;
end

% keeplines = 1:PFrate*dims(1); replacelines = PFrate*dims(1)+1:dims(1);
for iter=1:niters
    fprintf(' POCS iter %d\n',iter);
    % FT : convert image to kspace
    tmpim = fftshift(fft2(fftshift( (tmpim) )));
    % Replace acquired data for consistency
    tmpim(sampledPeLines,:) = kspace(sampledPeLines,:);
    % FT : convert kspace to image
    if useconj
    tmpim = conj(tmpim);
    tmpim = fftshift(fft2(fftshift(tmpim)));  
    else
    tmpim = prod(dims) * ifftshift(ifft2(ifftshift(tmpim)));  
    end
    tmpim = abs(tmpim) .* ph;

    if showprogress 
        tmp = (sqrt(sum(abs(tmpim(:,:,1,:).^2),4)));      % due to fftshift(), the 1st partition is the central one
        maxRange = sort( tmp(:), 'descend' );
        maxRange = maxRange( ceil(0.05 * numel(maxRange)) );    % ignore the "hottest" 5%
        if ~exist('pic','var')
            pic = [tmp tmp zeros(size(tmp),kprec)];
            diffScale = 1;
        else
            delta = abs( pic(:,dims(2)+(1:dims(2))) - tmp );
            diffScale = 0.5 * maxRange / median( delta(:) );
            pic(:,  dims(2)+(1:dims(2))) = tmp;
            pic(:,2*dims(2)+(1:dims(2))) = diffScale * delta;
            clear delta
        end
        imagesc(pic, [0 maxRange ]);
        title(sprintf('\\bfiteration %g\ninitial     |    current     |     abs(previous - current) × %g', iter, diffScale ))
        axis image
        colormap(gray(256))
        drawnow; 
        pause(plotpause);
        clear tmp
    end
end


%% POCS image reconstruction
pocs_image = tmpim;
pocs_kspace = fftshift(fft2(fftshift(tmpim)));

if makeplots

    f = figure('position', [163 86 429 780],'color','w');
    subplot(211);
    imshow(abs((pocs_image)),[]); title('PF+POCS (Magn)');
    set(gca,'FontSize',20); colorbar;
    subplot(212);
    % imshow(imag(pocs_recon_image),[]); title('Pocs recon (Phase)');
    imshow(angle(pocs_image),[]); title('PF+POCS (Phase)');
    set(gca,'FontSize',20,'Colormap',hsv); colorbar;
    
    
    
    
    diff_pf_pocs = image_fs - pocs_image;
    
    diff_pf_pocs_magn = abs(diff_pf_pocs);
    diff_pf_pocs_ph = angle(diff_pf_pocs);
    
    % diff_pf_pocs_magn = abs(image_fs)-abs(pocs_recon_image);
    % diff_pf_pocs_ph = angle(image_fs)-angle(pocs_recon_image);
    
    f = figure('position', [163 86 429 780],'color','w');
    subplot(211);
    imshow(diff_pf_pocs_magn,[]); title('Diff. b/w FS & PF(5/8)+POCS (Magn)');
    set(gca,'FontSize',20); colorbar;
    subplot(212);
    imshow(diff_pf_pocs_ph,[]); title('Diff. b/w FS & PF(5/8)+POCS (Phase)');
    set(gca,'FontSize',20,'Colormap',gray); colorbar;
    
    
    
    
    zp_image = ifftshift(ifft2(ifftshift( (kspace) )));
    diff_pf_zp = image_fs - zp_image;
    
    diff_pf_zp_magn = abs(diff_pf_zp);
    diff_pf_zp_ph = angle(diff_pf_zp);
    
    % diff_pf_pocs_magn = abs(image_fs)-abs(pocs_recon_image);
    % diff_pf_pocs_ph = angle(image_fs)-angle(pocs_recon_image);
    
    f = figure('position', [163 86 429 780],'color','w');
    subplot(211);
    imshow(diff_pf_zp_magn,[]); title('Diff. b/w FS & PF(5/8)+ZP (Magn)');
    set(gca,'FontSize',20); colorbar;
    subplot(212);
    imshow(diff_pf_zp_ph,[]); title('Diff. b/w FS & PF(5/8)+ZP (Phase)');
    set(gca,'FontSize',20,'Colormap',gray); colorbar;


    rmse.magn.zp = sqrt(mean( (abs(image_fs(:))-abs(zp_image(:))).^2 ) );
    rmse.magn.pocs = sqrt(mean( (abs(image_fs(:))-abs(pocs_image(:))).^2 ) );


end





