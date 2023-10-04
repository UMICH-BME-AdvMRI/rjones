

% Pick some values for T1, T2 (not specified in problem)
T1 = 300;	% ms.
T2 = 100;	% ms.

% Set specified flip angle
alpha = 60; % deg

% Pick some frequency range to eval ss response
df = 100; % Hz

ts = 1:1:1000;

nt = length(ts);

As0 = zeros(3,3,nt);
As = zeros(3,nt);
Bs = zeros(3,nt);

M0=[0 0 1]';

for ind=1:length(ts)
    t=ts(ind);
    [Afp,Bfp]=freeprecess(t,T1,T2,df);

    As(:,:,ind)=Afp*M0;
    Bs(:,ind)=Bfp;
end
