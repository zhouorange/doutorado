function out = prop_chan(prop)
%nargin=0;

if nargin==0,
 % parameters
 chan_type = 1;		% channel type
 Lam0 = 0.05e9;		% first cluster arrival rate (Hz)...set large to synchronize
 thresh_dB = -150;	% keep only paths whose powers exceed thresh_dB
 Lp = 100;		% number of paths per cluster
 % plotting
 plot_amp = 1;
 plot_pow = 1;

else%nargin>0
 chan_type = prop.chan_type;
 Lam0 = prop.Lam0;
 thresh_dB = prop.thresh_dB;
 Lp = prop.Lp;	
%Lbar = prop.Lbar;
%Lam = prop.Lam;
%lam1 = prop.lam1;
%lam2 = prop.lam2;
%beta = prop.beta;
%Gam = prop.Gam;
%sig_clust = prop.sig_clust;
%gam0 = prop.gam0;
%k_gam = prop.k_gam;
%m0 = prop.m0;
%m0hat = prop.m0hat;

 % plotting
 plot_amp = 0;
 plot_pow = 0;
end;%nargin


% set channel parameters based on type
if chan_type==1,
  Lbar = 13.6;		% avg number of clusters
  Lam = 0.0048e9; 	% cluster arrival rate (Hz)
  lam1 = 0.27e9; 	% inter-cluster arrival rate 1 (Hz)
  lam2 = 2.41e9; 	% inter-cluster arrival rate 2 (Hz)
  beta = 0.0078;	% inter-cluster arrival mixture coef
  Gam = 31.7e-9;	% cluster decay time
  gam0 = 3.7e-9;	% inter-cluster decay time at T=0
  k_gam = 0;		% inter-cluster decay dependence on T
  m0 = 0.77;		% for generating nakagami m-factor (dB)
  m0hat = 0.78;		% for generating nakagami m-factor (dB)
  sig_clust = 0;	% random power variation wrt exponential
elseif chan_type==2, 
  Lbar = 10.5;		% avg number of clusters
  Lam = 0.0243e9; 	% cluster arrival rate (Hz)
  lam1 = 0.15e9; 	% inter-cluster arrival rate 1 (Hz)
  lam2 = 1.13e9; 	% inter-cluster arrival rate 2 (Hz)
  beta = 0.062; 	% inter-cluster arrival mixture coef
  Gam = 104.7e-9;	% cluster decay time
  gam0 = 9.3e-9;	% inter-cluster decay time at T=0
  k_gam = 0;		% inter-cluster decay dependence on T
  m0 = 0.56;		% for generating nakagami m-factor (dB)
  m0hat = 0.25;		% for generating nakagami m-factor (dB)
  sig_clust = 0;	% random power variation wrt exponential
else
  error('chan_type not recognized')
end;

% generate number of clusters
L = max([1,poissrnd(Lbar)]);	% number of clusters

% generate cluster arrival times
T = cumsum([exprnd(1)/Lam0,exprnd(ones(1,L-1)/Lam)]);
gam = gam0 + k_gam*T;

% generate inter-cluster arrival times and powers
tau_ = zeros(L,Lp);
tau = [];
pow = [];
for l=1:L,
  dtau12 = [exprnd(ones(1,Lp-1)/lam1);exprnd(ones(1,Lp-1)/lam2)];
  tau_(l,:) = cumsum([0,dtau12(2-(rand(1,Lp-1)<beta))]);
  %tau_(l,:) = cumsum([0,exprnd(ones(1,Lp-1)/lam)]);
  tau = [tau, T(l)+tau_(l,:)];
  pow = [pow, exp(-T(l)/Gam - tau_(l,:)/gam(l) + sig_clust*randn(1,Lp)/10*log(10))/gam(l)/(lam1*(1-beta)+lam2*beta+1)];
end;
pow = pow/sum(pow);			% ensure unit-energy

% generate amplitudes
m = lognrnd(m0,m0hat)*ones(size(pow));	% lognormally distributed nakagami m-factors
amp = sqrt(gamrnd(m,pow./m));		% nakagami distributed amplitudes

% throw away negligible paths 
big = find(10*log10(pow)>thresh_dB);
Tau = tau(big);
Pow = pow(big);
Amp = amp(big);
Tau_max = T(1)-Gam*log(10^(thresh_dB/10));

% set outputs
out.tau = Tau;
out.amp = Amp;
out.pow = Pow;

% plot
if plot_pow,
 figure(1)
 plot(Tau,10*log10(Pow),'.');
 hold on; 
 plot(T,10*log10(exp(-(T-T(1))/Gam)),'g');
 for l=1:L,
   plot(T(l)+tau_(l,:),10*log10(exp(-(T(l)-T(1))/Gam-tau_(l,:)/gam(l))),'r');
 end;
 hold off;
 grid on;
 axis([0,Tau_max,thresh_dB,0]);
 ylabel('power [dB]');
 xlabel('delay [sec]');
end;%plot_pow

if plot_amp,
 figure(2)
 subplot(211)
   handy = stem(Tau,Amp); set(handy,'MarkerSize',0);
   grid on
   ylabel('amplitude');
 subplot(212)
   handy = stem(Tau,20*log10(abs(Amp))); set(handy,'MarkerSize',0);
   grid on
   ylabel('amplitude [dB]');
 xlabel('delay [sec]');
end;%plot_amp
