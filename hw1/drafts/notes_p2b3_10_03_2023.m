%% vv
% Bloch Equation Simulation, Excercise B-5b
% -----------------------------------------
% 

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 2;		% ms.
TR = 10;
flip = [0:0.01:0.5 ]*pi;	% radians. % [0:0.01:0.5 ]*pi;	% radians.
inc = 117/180*pi;
Nex = 100;

Ms = zeros(length(flip),3);
for k=1:length(flip)
	sig1(k)=spgrsignal(flip(k),T1,T2,TE,TR,df,Nex,inc,100,1,1);
    [sig1,Mss,Msigs,Ms]=spgrsignal(flip(k),T1,T2,TE,TR,df,Nex,inc,100,1,0);
    Ms(k,:)=Mss;

	[Msig,M]=srsignal(flip(k),T1,T2,TE,TR,df);
	sig2(k)=M(1)+i*M(2);
end;


figure; plot(Ms);

% ===== Plot the Results ======

plot(flip,abs(sig1),'r-',flip,abs(sig2),'b--');
xlabel('Flip (rad)');
ylabel('Signal Magnitude');
grid on;
title('Signal vs Flip for RF-spoiled GRE and SR');
legend('RF-spoiled GRE','Saturation-Recovery Approximation');



%% ww
% % Bloch Equation Simulation, Excercise B-1d
% % -----------------------------------------
% % 
% 
% df = 0;		% Hz off-resonance.
% T1 = 600;	% ms.
% T2 = 100;	% ms.
% TE = 1;		% ms.
% TR = 500;	% ms.
% flip = pi/3;	% radians.
% 
% 
% M = [0;0;1];
% Rflip = yrot(flip);
% [Atr,Btr] = freeprecess(TR,T1,T2,df);
% 
% % M1 = Atr * Rflip * M + Btr.
% %
% % But M1=M in steady state, so
% %
% % 	M = Atr*Rflip * M + Btr.
% %	(I-Atr*Rflip)*M = Btr.
% 
% M = inv(eye(3)-Atr*Rflip)*Btr