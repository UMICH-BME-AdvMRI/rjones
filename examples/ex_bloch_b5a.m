% Bloch Equation Simulation, Excercise B-5a
% -----------------------------------------
% 

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 2;		% ms.
TR = 10;
flip = pi/6;	% radians.
inc = 117/180*pi;

Nex = 100;

Nf = 100;	% Simulate 50 different gradient-spoiled spins.
phi = [1:Nf]/Nf*2*pi;

M=zeros(3,Nf,Nex+1);
Msig = zeros(Nex,1);


	
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

	for k=1:Nf
		M(:,k) = zrot(phi(k))*M(:,k);
	end;

	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end;


% ===== Plot the Results ======

time = [0:Nex-1]*TR+TE;
subplot(2,1,1);
plot(time,abs(Msig));
xlabel('Time (ms)');
ylabel('Magnitude');
grid on;


subplot(2,1,2);
plot(time,angle(Msig));
xlabel('Time (ms)');
ylabel('Phase');
grid on;

		
		






