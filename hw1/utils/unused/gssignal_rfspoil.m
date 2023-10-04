% 
%	function [Msig,Mss] = gssignal(flip,T1,T2,TE,TR,dfreq,phi)
% 
%	Calculate the steady state signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  dfreq is the resonant
%	frequency in Hz.  flip is in radians.
%	phi is the phase twist at the end of the sequence.

function [Msig,Mss] = gssignal_rfspoil(flip,T1,T2,TE,TR,dfreq,phi,theta)

Rflip = throt(flip,theta);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);
[Ate,Bte] = freeprecess(TE,T1,T2,dfreq);

% 	To add the gradient spoiler twist, we just
%	multiply Atr by zrot(phi):

Atr = zrot(phi)*Atr;


%%%%
[Ate,Bte] = freeprecess(TE,T1,T2,0);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,0);
M = [zeros(2,Nf);ones(1,Nf)];
on = ones(1,Nf);	
Rfph = 0;
Rfinc = inc;
for n=1:Nex
	A = Ate * throt(flip,Rfph);
	B = Bte;
	M = A*M+B*on;
	Msig(n) = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*Rfph);	
	M=Atr*M+Btr*on;
	for k=1:Nf, M(:,k) = zrot(phi(k))*M(:,k); end;
	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end;
%%%%

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

Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
Msig = Mss(1)+i*Mss(2);

