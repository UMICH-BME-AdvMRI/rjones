function [Msig,Mss] = gssignal(flip,T1,T2,TE,TR,dfreq,phi)
%	function [Msig,Mss] = gssignal(flip,T1,T2,TE,TR,dfreq,phi)
% 
%	Calculate the steady state signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  dfreq is the resonant
%	frequency in Hz.  flip is in radians.
%	phi is the phase twist at the end of the sequence.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Author: Brian Hargreaves 
%
% ====REFERENCE==== ! PLEASE NOTE !
%  This function was from Brian Hargreaves' Bloch Equation Simulation webpage:
%   http://mrsrl.stanford.edu/~brian/bloch/


%  Apply RF excitation - tip spins into transverse plane
Rflip = yrot(flip);

%  Simulate spin precession w Bloch eqns
[Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);
[Ate,Bte] = freeprecess(TE,T1,T2,dfreq);

% 	To add the gradient spoiler twist, we just
%	multiply Atr by zrot(phi):
Atr = zrot(phi)*Atr;

% Let 	M1 be the magnetization just before the tip.
%	M2 be just after the tip.
%	M3 be at TE.
%
% then
%	M2 = Rflip * M1
%	M3 = Ate * M2 + Bte
%	M1 = Atr * M3 + Btr
%
% Solve for M3...
%
%	M3 = Ate*Rflip*Atr*M3 + (Ate*Rflip*Btr+Bte)

% Calculate magnetization & steady state signal at TE
Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
Msig = Mss(1)+1i*Mss(2);

