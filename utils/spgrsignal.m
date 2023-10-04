%	function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc)
%
%	Function calculates the signal from an RF-spoiled sequence
%	after Nex excitations.
%
function [Msig,Mss,Msigs]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc,Nf,gs_cycles)

if nargin<10
    gs_cycles = 1;
end
if nargin<9
    Nf = 100;
end
if (nargin < 8)
	inc = 117/180*pi;
end
if (nargin < 7)
	Nex = 100;
end
if (nargin < 6)
	df = 0;
end

% Nf = 100;	% Simulate 100 different gradient-spoiled spins.
phi = [1:Nf]/Nf*2*pi*gs_cycles;

M=zeros(3,Nf,Nex+1);
Msigs = zeros(Nex,1);
	
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);

M = [zeros(2,Nf);ones(1,Nf)];
on = ones(1,Nf);
	
Rfph = 0;
Rfinc = inc;

for n=1:Nex

	A = Ate * throt(flip,Rfph);
	B = Bte;
	M = A*M+B*on;

	Msig = mean( squeeze(M(1,:)+1i*M(2,:)) ) * exp(-1i*Rfph);
    Msigs(n)=Msig;
	Mss = M;

	M=Atr*M+Btr*on;

	for k=1:Nf
		M(:,k) = zrot(phi(k))*M(:,k);
    end

	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end


		






