function [Msig,Mss,Msigs,Ms]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc,Nf,gs_cycles,MakePlots)
%	function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc,Nf,gs_cycles)
%
%	Function calculates the signal from an RF-spoiled sequence
%	after Nex excitations.
% INPUTS:
%   flip        = flip angle (rad)
%   T1,T2,TE,TR = ms   
%   df          = off-resonance frequency (Hz)
%   Nex         = # of excitations (RF pulses) to apply for each simulation
%   inc         = increment to add to RF pulse phase after each excitation
%   Nf          = # of gradient spoiling frequencies to simulate in each
%                   voxel
%   gs_cycles   = # of gradient spoiling dephasing cycles (2*pi rad, e.g.,
%                   1cycle=2pi dephasing)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Author: Brian Hargreaves 
%
% PLEASE NOTE
%   This function was from Brian Hargreaves Bloch Equation Simulation webpage
%   http://mrsrl.stanford.edu/~brian/bloch/


if nargin<11
    MakePlots=0;
end
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

Ms=zeros(3,Nf,Nex+1);
Msigs = zeros(Nex,1);
MsigsExc = zeros(Nf,Nex);
	
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

    Ms(:,:,n) = M;
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


		
if MakePlots
    % ===== Plot the Results ======

    f=figure;
    f.Name = 'spgrsignal';

    time = [0:Nex-1]*TR+TE;
    subplot(2,1,1);
    plot(time,abs(Msigs));
    xlabel('Time (ms)');
    ylabel('Magnitude');
    grid on;
    
    
    subplot(2,1,2);
    plot(time,angle(Msigs));
    xlabel('Time (ms)');
    ylabel('Phase');
    grid on;

    drawnow;

end





