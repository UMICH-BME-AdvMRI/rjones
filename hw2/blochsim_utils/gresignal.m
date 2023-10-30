function [Msig,Mss] = gresignal(flip,T1,T2,TE,TR,dfreq,ncycles,N)
%	function [Msig,Mss] = gresignal(flip,T1,T2,TE,TR,dfreq)
% 
%	-Calculate the steady state gradient-spoiled signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  
%  -dfreq is the resonant frequency in Hz.  flip is in radians.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Author: Brian Hargreaves 
%
% PLEASE NOTE
%   This function was from Brian Hargreaves Bloch Equation Simulation webpage
%   http://mrsrl.stanford.edu/~brian/bloch/




% N = 100;
M = zeros(3,N);
% phi = linspace(-0.5,0.5,N) * ncycles*pi; %([1:N]/N-0.5 )
phi = ([1:N]/N-0.5 ) * ncycles*2*pi;

for k=1:N
	[M1sig,M1] = gssignal(flip,T1,T2,TE,TR,dfreq,phi(k));
	M(:,k)=M1;
end


Mss = mean(M')';
Msig = Mss(1)+i*Mss(2);



