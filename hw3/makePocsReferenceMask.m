function [ mask ] = makePocsReferenceMask( dims, pf )
% [ mask ] = makePocsReferenceMask( dims, pf )
% 
% create kspace sampling mask for pocs reference (low freq, fully sampled)
%
% inputs:
%   dims    = dimensions of kx-ky (2-element vector)
%   pf      = partial fourier factor (fraction/float in interval (0,1)
% outputs:
%   mask    = mask with 1s where pocs reference is used and 0s where its not

if nargin<2
    error('need 2 input args');
end
assert(length(dims)==2);

DCoffset = ceil(dims(1)/2)+1;
dcradius = round(dims(1)*(pf-0.5));
mask = zeros(dims);
mask(DCoffset-dcradius:DCoffset+dcradius-1,:)=1;
