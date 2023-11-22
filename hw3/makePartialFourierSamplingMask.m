function [ mask ] = makePartialFourierSamplingMask( dims, pf )
% [ mask ] = makePartialFourierSamplingMask( dims, pf )
% 
% create kspace sampling mask for partial fourier acq
%
% inputs:
%   dims    = dimensions of kx-ky (2-element vector)
%   pf      = partial fourier factor (fraction/float in interval (0,1)
% outputs:
%   mask    = mask with 1s where kspace is acquired and 0s where its not

if nargin<2
    error('need 2 input args');
end
assert(length(dims)==2);
mask = zeros(dims);
mask(1:round(pf*dims(1)),:)=1;

