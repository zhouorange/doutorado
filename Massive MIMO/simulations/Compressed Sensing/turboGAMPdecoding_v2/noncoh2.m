function out = noncoh2(run,ofdm,chan,ldpc,eq)
%nargin=0;

if nargin==0,

 % OFDM and channel params
 N = 1024;		% # subcarriers
 Lx = round(N/4);	% impulse response length
M = 4;			% # bits per subcarrier (must be even)
SNRdB = 14.0;		% SNR in dB (note: SNR = bpcu * Eb/No)
Np = round(0.875*4*64);	% # pilot-only subcarriers [Nyquist=256=4*64, dflt=224=0.875*4*64]
Mt = round(0*64*M);	% # training coded-bits (per OFDM symbol) [dflt=1.75*64*M]
 uwb = 1; 		% use UWB channel model?
 uwb_type = 2;		% UWB: 1=outdoor-LOS, 2=outdoor-NLOS [dflt=2]
 alf_rc = 0.5;		% UWB: raised-cosine parameter [dflt=0.5]
 baud = 1.0e-6/Lx;	% UWB: baud interval [dflt=1.0e-6/Lx]
 Lam0 = 1/baud;		% UWB: first-arrival rate (Hz) [dflt=1/baud ]
 Lpre = 20;		% UWB: number of precursor taps [dflt=20]
 pi01 = 0.025;		% nonUWB: support switching prob [dflt=0.1]
 pi10 = 0.00625;	% nonUWB: support switching prob [dflt=0.025]
 clust = 10;		% nonUWB: initial cluster length
 Le = 15;		% nonUWB: half-power length for exponential PDP [dflt=13]
 rho = 1e-3;		% nonUWB: off/on tap-variance ratio (0<=rho<1) [dflt=1e-5]

 % LDPC params
 Mc_target = 8192;	% target # coded bits (total) [dflt=8192]
bpcu_target = 2.0;	% target bits-per-channel-use (including pilots & training bits)
 weight = 3;		% mean column weight [dflt=3]
 seed = 1;		% code generation seed (MATLAB 5 uniform generator) [dflt=1]
 LLRmax_ldpc = 20;	% LDPC decoder max posterior LLR [dflt=20]
 iter_ldpc = 25;	% LDPC decoding iterations [dflt=25]
 interleave = 1;	% use random interleaving? [dflt=1]

 % turbo parameters
 stop_ber = 1;		% 1:stop when ber=0, 2:stop when ldpc successful
 LLR_diff_stop = 1e-2;	% stops turbo below this average LLR difference [dflt=1e-1 or 1e-2?]
 LLRmax_eq = 10;	% equalizer max extrinisic bit-LLR [dflt=10]
 iter_turbo = 10;	% max turbo iterations (>0) [dflt=10]

 % equalizer parameters
 iid_bg = 0; 		% use iid BG model instead of nid GM2 model [dflt=0]
 iter_mc = 5;		% markov-chain iterations (per turbo iteration) [dflt=0,5] 
 iter_amp = 15; 	% max GAMP iterations [dflt=15]
 llrmax_eq = 10;	% equalizer max extrinisic support-llr [dflt=10]
 llrmax_mc = 10;	% markov-chain decoder max extrinisic support-llr [dflt=10]
 ratio_ez = 0.99;	% max allowed mu_e/mu_z in GAMP (< 1) [dflt=0.99]
 reg_amp = 0.02;	% GAMP regularization factor; helps stability at high SNR [dflt=0.02]
 llr_diff_stop = 1e-2;	% stops MC iters below this average llr difference [dflt=1e-1 or 1e-2?]
 amp_diff_stop = 1e-6;	% stops GAMP below this average energy difference [dflt=1e-6]
 %spgl1_scale = [1.5:-0.2:0.7];     % list of SPGL1 thresholds [dflt=[1.1:-0.1:0.7]]
 spgl1_scale = [1.5,0.9];     % list of SPGL1 thresholds [dflt=[1.5,0.9]]
 spgl1_opts = spgSetParms('verbosity',0,'iterations',N); % SPGL1 options
 hard_fdbk = 0;		% use hard feedback for turbo LMMSE and turbo LASSO? [dflt=0]

 % plotting params
 plot_alph = 0;		% figure number of plot (or zero for no plot)
 plot_chan = 1;		% figure number of plot (or zero for no plot)
 plot_amp_eq = 2;	% figure number of plot (or zero for no plot)
 plot_amp_dec = 3;	% figure number of plot (or zero for no plot)
 plot_lin_eq = 4;	% figure number of plot (or zero for no plot)
 plot_lin_dec = 5;	% figure number of plot (or zero for no plot)
 plot_las_eq = 6;	% figure number of plot (or zero for no plot)
 plot_las_dec = 7;	% figure number of plot (or zero for no plot)

 % run parameters
turbo_amp = 1;		% run turbo GAMP 
genie_est = 0;		% run genie channel estimation w/ coherent decoding
turbo_lin = 0;		% run turbo LMMSE
turbo_las = 0;		% run turbo LASSO
 oneshot_lasso = 0;	% run one-shot LASSO
 genie_debug = 0; 	% uses perfectly known info bits (only used for debugging!)

 % set randn and rand seeds
 %reset(RandStream.getDefaultStream);	% uncomment to start at reference
 %defaultStream.State = savedState;	% uncomment to use previous state
 defaultStream = RandStream.getDefaultStream;	% leave uncommented!
 savedState = defaultStream.State;	% leave uncommented!

elseif nargin==5,

 genie_est = run.genie; 
 turbo_amp = run.amp;
 turbo_lin = run.lin;
 turbo_las = run.las;
 oneshot_lasso = 0;
 N = ofdm.subcar;
 Np = ofdm.pilot_subcar;
 M = ofdm.bits_per_qamsym;
 Mt = ofdm.trainingbits_per_ofdmsym;
 Lx = chan.len;
 Le = chan.halfpowerlen;	% non-uwb
 clust = chan.clust;		% non-uwb
 pi10 = chan.pi10;		% non-uwb
 pi01 = chan.pi01;		% non-uwb 
 rho = chan.onoffratio;		% non-uwb
 uwb = chan.uwb;		% uwb
 uwb_type = chan.uwb_type;	% uwb
 alf_rc = chan.alf_rc;		% uwb
 baud = chan.baud;		% uwb
 Lam0 = chan.Lam0;		% uwb
 Lpre = chan.Lpre;		% uwb
 SNRdB = chan.SNRdB;
 Mc_target = ldpc.codelen_target;
 bpcu_target = ldpc.bpcu_target;
 weight = ldpc.weight;
 seed = ldpc.seed;
 LLRmax_ldpc = ldpc.LLRmax;
 iter_ldpc = ldpc.iter;
 interleave = ldpc.interleave;
 LLRmax_eq = eq.LLRmax;
 iid_bg = eq.iid_bg;
 if turbo_amp||turbo_lin||turbo_las,
   iter_turbo = eq.iter_turbo;
   stop_ber = eq.stop_ber;
   LLR_diff_stop = eq.LLR_diff_stop;
 end;
 if turbo_amp,
   reg_amp = eq.reg_amp;
   ratio_ez = eq.ratio_ez;
   iter_mc = eq.iter_mc;
   llrmax_eq = eq.llrmax_eq;
   llrmax_mc = eq.llrmax_mc;
   iter_amp = eq.iter_amp;
   llr_diff_stop = eq.llr_diff_stop;
   amp_diff_stop = eq.amp_diff_stop;
 end;
 if turbo_las||oneshot_lasso,
   spgl1_scale = eq.spgl1_scale;
   spgl1_opts = eq.spgl1_opts;
 end;
 plot_alph = 0;
 plot_chan = 0;
 plot_amp_eq = 0;
 plot_amp_dec = 0;
 plot_lin_eq = 0;
 plot_lin_dec = 0;
 plot_las_eq = 0;
 plot_las_dec = 0;
 genie_debug = 0; 

else
  error('incorrect number of parameters')
end;

 % calculate channel quantities
 sig2v = 10^(-SNRdB/10);		% freq-domain AWGN variance
 mu_v = (sig2v+reg_amp)*ones(N,1);	% apriori noise vars (need to add 0.01 or 0.02!)
 if uwb,
   if iid_bg, bggm_str = 'bgstats_type='; else bggm_str = 'gmstats_type='; end;
   file_str = [bggm_str,num2str(uwb_type),...
   		'_L=',num2str(Lx),...
		'_TL=',num2str(baud*Lx*1e6),...
		'_alf=',num2str(alf_rc),...
		'_TLam=',num2str(baud*Lam0),'.mat'];
   if exist(file_str)==2,
     load(file_str);			% provides lam,mu_on,mu_off,p01,p10,E_h
   else
     error(['You must first configure and run ofdm_chan.m to generate the file ',file_str])
   end;
   prop.chan_type = uwb_type;  		% uwb type
   prop.Lam0 = Lam0;      		% first cluster arrival rate (Hz)
   prop.thresh_dB = -150; 		% model precision [dflt=-150]
   prop.Lp = 100;         		% paths per cluster [dflt=100]
 else
   p01 = [zeros(1,clust)+eps,pi01*ones(1,Lx-clust-1)];% sparsity pattern switching prob
   p10 = [ones(1,clust)-eps,pi10*ones(1,Lx-clust-1)];% sparsity pattern switching prob
   lam = [1;NaN*ones(Lx-1,1)]; 
   lam_ = [1,0];
   for l=1:Lx-1,
     Pi=[1-p01(l),p01(l);p10(l),1-p10(l)];	% state transition mtx 
     lam_ = lam_*Pi;  
     lam(l+1) = lam_(1);
   end;
%  p01 = pi01*ones(1,Lx-1);		% sparsity pattern switching prob 
%  p10 = pi10*ones(1,Lx-1);		% sparsity pattern switching prob
   %if Ls<1, error('Ls<1!'); end;	% need at least one sparse channel tap
   %lam = Ls/Lx*ones(Lx,1);		% channel activity factor
   dpp = 2.^(-[0:Lx-1]/Le);  dpp = (dpp.')/mean(dpp); % normalized PDP
   mu_on = (dpp/Lx)./(lam+(1-lam)*rho);	% apriori on-tap var 
   mu_off = rho*mu_on;			% apriori off-tap var
 end;
 mu_onoff = lam.*mu_on+(1-lam).*mu_off;	% apriori tap var (sums to 1) <=> PDP
 gam = p10(:)./lam(1:Lx-1); 		% markov independence parameter

 % calculate modulation quantities
 if M~=2*floor(M/2), error('M must be even!'); end;
 R_target = N/((N-Np)*M-Mt)*bpcu_target;
 if R_target>=1, error('need R_target<1!'); end;
 Phi = dftmtx(N); Phi = Phi(:,1:Lx);	% NxN DFT matrix (not unitary!)
 Nd = N-Np;				% # data/training subcarriers
 T = round(Mc_target/(Nd*M-Mt));	% # OFDM symbols (total)
 Mb = T*Nd*M;				% # coded data/training bits (total)
 Mc = Mb-T*Mt;				% # coded data bits (total)
 Mi = round(Mc*R_target);		% # info bits (total)
 R = Mi/Mc;
 bpcu = R*((N-Np)*M-Mt)/N;
 EbNo = 1/bpcu/sig2v;

 % eventually export as outputs
 clear out;
 out.Mc = Mc;
 out.R = R;
 out.T = T;

 % generate channel and noise realizations
 x_true = cell(1,T);				% impulse response
 llr_true = cell(1,T);				% support llrs
 supp_true = cell(1,T);				% support pattern
 comp_true = cell(1,T);				% support pattern complement
 mu_true = cell(1,T);				% per-tap variance
 z_true = cell(1,T);				% freq response
 v_true = cell(1,T);				% noise samples 
 for t=1:T,
  if uwb,					% use UWB realizations
   mpc = prop_chan(prop);	% generate uwb delays & amplitudes
   %mpc.tau = i/M*T; mpc.amp = 1;% trivial channels for testing
   n_s = [0:Lx-1]'+eps;		% baud-normalized sampling times
   x_true{t} = zeros(Lx,1);	% sampled channel
   for p=1:length(mpc.tau),	% for each path... 
     if mpc.tau(p)/baud>Lx, break; end;
     n_p = Lpre+mpc.tau(p)/baud;% baud-normalized propagation delay
     x_true{t} = x_true{t} + exp(1i*2*pi*rand(1))*mpc.amp(p)*cos(...
        pi*alf_rc*(n_s-n_p))./(1-(2*alf_rc*(n_s-n_p)).^2).*sinc(n_s-n_p);
   end;
   x_true{t} = x_true{t}/sqrt(E_h);		% normalize to unit energy
   %CNoff = exp(-(abs(x_true{t}).^2)./mu_off)./mu_off;
   %CNon = exp(-(abs(x_true{t}).^2)./mu_on)./mu_on;
   %gamon = (lam.*CNon + eps)./((1-lam).*CNoff+lam.*CNon + eps); %posterior on supp
   llr_true{t} = (abs(x_true{t}).^2).*(mu_on-mu_off)./mu_on./mu_off -log(mu_on)+log(mu_off);
   supp_true{t} = find(llr_true{t} - log(1./lam-1) > 0);
   comp_true{t} = setdiff([1:Lx],supp_true{t});	
   mu_true{t} = zeros(Lx,1);
   mu_true{t}(supp_true{t}) = mu_on(supp_true{t});
   mu_true{t}(comp_true{t}) = mu_off(comp_true{t});
  else						% use synthetic BG realizations 
   supp_true{t} = [];
   while abs(mean(lam)-length(supp_true{t})/Lx)/mean(lam) > 0.2, % make sure lam is matched!
     pattern = NaN*ones(Lx,1);			% binary support pattern
     pattern(1) = (rand(1)<lam(1));		% initial state
     for l=1:Lx-1,				% run Markov chain
       pattern(l+1) = (rand(1) < [pattern(l),1-pattern(l)]*[1-p01(l);p10(l)]);
     end;
     supp_true{t} = find(pattern);
   end;
   comp_true{t} = setdiff([1:Lx],supp_true{t});
   mu_true{t} = zeros(Lx,1);
   mu_true{t}(supp_true{t}) = mu_on(supp_true{t});
   mu_true{t}(comp_true{t}) = mu_off(comp_true{t});
   x_true{t} = zeros(Lx,1);
   x_true{t} = sqrt(mu_true{t}).*(randn(Lx,2)*[1;1i]/sqrt(2));
  end;
  z_true{t} = fft(x_true{t},N);
  v_true{t} = randn(N,2)*[1;1i]*sqrt(sig2v/2);	% noise samples
%v_true{t} = zeros(N,1);

  if plot_chan&&(t==1),
    figure(plot_chan); 
    subplot(111)
    handy = plot(20*log10(abs(x_true{t})),'g');
    hold on; 
     handy = [plot(supp_true{t},20*log10(abs(x_true{t}(supp_true{t}))),'b.');handy];
     handy = [handy;plot(10*log10(mu_on),'b-.')];
     handy = [handy;plot(10*log10(mu_onoff),'k')];
     handy = [handy;plot(10*log10(max(0,log((mu_off+eps)./mu_on./(1./lam -1)...
					)./(1./mu_on-1./(mu_off+eps)))),'k:')];
     handy = [handy;plot(10*log10(mu_off+eps),'r--')];
     handy = [handy;plot(comp_true{t},20*log10(abs(x_true{t}(comp_true{t}))),'ro')]; 
    hold off;
    set(handy(1),'Markersize',13)
    set(handy(7),'Markersize',4)
    ylabel('dB')
    xlabel('lag [baud]')
    title('channel realization')
    axe=axis; axis([0,Lx,axe(3:4)])
    legend(handy,'taps: big','channel','var: big','PDP','threshold','var: small','taps: small','Location','Best')
    drawnow;
  end;%plot_chan
 end;
%for t=1:T, E_x = norm(x_true{t})^2, end;		% uncomment to check energy
%for t=1:T, lam_x = length(supp_true{t})/Lx, end;	% uncomment to check sparsity

 % build pilot/training pattern
 if Np==0,
   n_pilot_ = [];
 else
   %n_pilot_ = floor([1:N/Np:N]);	% pilot subcarrier indices, in {1...N}
   %n_pilot_ = floor([1+floor(N/Np/2):N/Np:N]);	% pilot-subcarrier indices, in {1...N}
   tmp = randperm(N); n_pilot_ = sort(tmp(1:Np));
 end;
 %k_pilot_ = ones(Np,1);		% not fair: raises avg signal power!
 k_pilot_ = ceil(rand(Np,1)*2^M);	% pilot alphabet indices, in {1...2^M}
 n_data_ = setdiff([1:N],n_pilot_);	% data/training subcarrier indices, in {1...N}
 if Mt>Nd*M, error('too many training bits!'); end;
 if Mt==0,				% if no training bits...
   m_train_ = [];
 elseif Mt<=2*Nd			% if all training bits can be MSBs...
   Nmsb = floor(Mt/2);			% # subcarriers whose MSBs are training bits
   m_train_ = sort([M*floor(0:Nd/Nmsb:Nd-1)+1, M*floor(0:Nd/Nmsb:Nd-1)+2]);
 else
   m_train_ = round([1:Nd*M/Mt:Nd*M]);	% indices of training bits, in {1...Nd*M}
   display('warning: training bits are not MSBs!')
 end;
 mm_train_ = []; 			% indices of training bits, in {1...Mb}
 for t=1:T,
   mm_train_ = [mm_train_,(t-1)*Nd*M+m_train_]; 
 end;
 mm_data_ = setdiff([1:Mb],mm_train_);	% indices of data bits, in {1...Mb}
 bb_train  = ones(1,T*Mt);		% values of training bits (better than random!)
 %bb_train  = round(rand(1,T*Mt));	% values of training bits

 % generate symbol and bit hypotheses, assuming multi-level Gray mapping
 B_hyp = [0;1];
 for m=1:M-1,
   B_hyp=[zeros(2^m,1),B_hyp;ones(2^m,1),B_hyp]; % bit-seq hypos 
 end;
 zB_hyp = [zeros(2^M,2),B_hyp];
 s_hyp = zeros(2^M,1);			% symbol alphabet
 for q=1:M/2,
   new = (2*zB_hyp(:,[2*q+1,2*q+2])-1)*[1;-1i]/2^q;
   for l=1:q,
     new = ((-1).^zB_hyp(:,2*l-1)).*real(new) ...
                 + 1i*((-1).^zB_hyp(:,2*l)).*imag(new);
   end;
   s_hyp = s_hyp + new;
 end;
 s_hyp = (s_hyp/sqrt(mean(abs(s_hyp).^2))).';	
 s_hyp_inv = 1./s_hyp;
 s_hyp2 = abs(s_hyp).^2;
 k_one = zeros(M,2^(M-1)); 		% alphabet indices with m'th bit=1
 for m=1:M, k_one(m,:) = find(B_hyp(:,m)==1).'; end;
 k_zero = zeros(M,2^(M-1)); 		% alphabet indices with m'th bit=1
 for m=1:M, k_zero(m,:) = find(B_hyp(:,m)==0).'; end;
 if plot_alph,
   figure(plot_alph);
   subplot(111)
   plot(real(s_hyp),imag(s_hyp),'.')
   for k=1:2^M,
     handy = text(real(s_hyp(k))-0.04,imag(s_hyp(k))-0.03,num2str(B_hyp(k,:)));
     set(handy,'FontSize',6);
   end;
 end;% plot_alph

 % load/generate LDPC matrices
 str_ldpc = ['./pargen/pargen_',num2str(Mc-Mi),'_',num2str(Mc),...
 	'_',num2str(weight),'_',num2str(seed),'.mat'];
 if exist(str_ldpc,'file'),
   load(str_ldpc);
 else
%  if strfind(computer('arch'),'64'),
%    error('Unable to generate these LDPC matrices with 64-bit Matlab!');
%  else
     H = ldpc_generate1(Mc-Mi,Mc,weight,2,seed);
     [HH, G] = ldpc_h2g1(H);			% systematic: G = [A|I]
     save(str_ldpc,'HH','G');
%  end;
 end;

 % generate coded bits and interleave them
 ii_true = (rand(Mi,1) > 0.5);			% info bits 
 cc_true = mod(ii_true'*G,2)';			% coded bits 
 if interleave,
   intlv = randperm(Mc);			% interleave 
 else,
   intlv = [1:Mc]; 
 end;
 bb_true = NaN*ones(Mb,1);			% interleaved/coded + training bits
 bb_true(mm_train_) = bb_train;
 bb_true(mm_data_) = cc_true(intlv);
 indx_bb1 = find(bb_true==1);
 indx_bb0 = find(bb_true==0);

 % generate Gray-mapped symbols from interleaved&training bits, and add pilot symbols
 s_true = cell(1,T);
 bin_2_int = 2.^[M-1:-1:0];
 for t=0:T-1,
   s_true{t+1} = NaN*ones(N,1);
   s_true{t+1}(n_pilot_) = s_hyp(k_pilot_).'; 
   for i=0:Nd-1,
     s_true{t+1}(n_data_(i+1)) = s_hyp(1+bin_2_int*bb_true(t*Nd*M+i*M+[1:M])).';
   end;
 end;

 % build noisy OFDM observations
 y = cell(1,T);
 for t=1:T,
   y{t} = z_true{t}.*s_true{t} + v_true{t};
 end;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Genie Channel Estimation %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if genie_est,					% try various genie chan estimators

  % generate observations with different noise realization for bit-aware genie 
  v_fake = cell(1,T);
  y_fake = cell(1,T);
  for t=1:T,
    v_fake{t} = randn(N,2)*[1;1i]*sqrt(sig2v/2);
    y_fake{t} = z_true{t}.*s_true{t} + v_fake{t};
  end;

  % bits -> symbols
  p_s_to_y = cell(1,T);			
  LLR_bb_to_del = zeros(Mb,1);			% initialize bit LLRs 
  LLR_bb_to_del(mm_train_) = -log(eps)*(2*bb_train-1); % training bits
  for t=1:T,
    p_b_to_del = 1./(1+exp(-LLR_bb_to_del((t-1)*M*Nd+[1:Nd*M]))); % LLRs -> bit probs
    p_s_to_y{t} = zeros(N,2^M);
    for i=0:Np-1,				% pilot symbol probs
      p_s_to_y{t}(n_pilot_(i+1),k_pilot_(i+1)) = 1;	 
    end;
    for i=0:Nd-1,				% data/training symbol probs
      p_Bi = ones(2^M,1)*p_b_to_del(i*M+[1:M]).';
      p_s_to_y{t}(n_data_(i+1),:) = prod(B_hyp.*p_Bi+(1-B_hyp).*(1-p_Bi),2).';
    end;
  end;%t

  % channel estimation
  nmse_bsg = NaN*ones(T,1);
  nmse_bgl = NaN*ones(T,1);
  nmse_sg = NaN*ones(T,1);
  x_hat_bsg = cell(1,T); z_hat_bsg = cell(1,T); mu_z_bsg = cell(1,T);
  x_hat_bgl = cell(1,T); z_hat_bgl = cell(1,T); mu_z_bgl = cell(1,T);
  x_hat_sg = cell(1,T); z_hat_sg = cell(1,T); mu_z_sg = cell(1,T);
  mu_z_pcsi = cell(1,T);
  for t=1:T,
    A = diag(s_true{t})*Phi;
    % perfect-CSI genie
    mu_z_pcsi{t} = eps*ones(N,1);
    % bit & support aware genie (non-compressed sensing)
    %x_hat_bsg{t} = mu_true{t}.*(A'*((A*diag(mu_true{t})*A'+sig2v*eye(N))\y_fake{t}));
    x_hat_bsg{t} = mu_true{t}.*(A'*((A*diag(mu_true{t})*A'+sig2v*eye(N))\y{t}));
    z_hat_bsg{t} = fft(x_hat_bsg{t},N);
    mu_z_bsg{t} = ones(N,1)*(norm(z_hat_bsg{t}-z_true{t})^2)/N;
    nmse_bsg(t) = norm(x_hat_bsg{t}-x_true{t})^2/norm(x_true{t})^2;
    % bit-aware genie LMMSE (non-compressed sensing)
    x_hat_bgl{t} = mu_onoff.*(A'*((A*diag(mu_onoff)*A'+sig2v*eye(N))\y_fake{t}));
    %x_hat_bgl{t} = mu_onoff.*(A'*((A*diag(mu_onoff)*A'+sig2v*eye(N))\y{t}));
    z_hat_bgl{t} = fft(x_hat_bgl{t},N);
    mu_z_bgl{t} = ones(N,1)*(norm(z_hat_bgl{t}-z_true{t})^2)/N;
    nmse_bgl(t) = norm(x_hat_bgl{t}-x_true{t})^2/norm(x_true{t})^2;
    % compressed channel sensing
    if Np>0, 
      A = diag(s_true{t}(n_pilot_))*Phi(n_pilot_,:);		
      % support-aware genie (i.e., pilot-aided)
      x_hat_sg{t} = mu_true{t}.*( ...
      	A'*((A*diag(mu_true{t})*A'+sig2v*eye(Np))\y{t}(n_pilot_)));
      z_hat_sg{t} = fft(x_hat_sg{t},N);
      mu_z_sg{t} = ones(N,1)*(norm(z_hat_sg{t}-z_true{t})^2)/N;
      nmse_sg(t) = norm(x_hat_sg{t}-x_true{t})^2/norm(x_true{t})^2;
    end;%if Np>0
  end;%t=1:T
  nmse_bsg_dB = 10*log10(mean(nmse_bsg));
  nmse_bgl_dB = 10*log10(mean(nmse_bgl));
  nmse_sg_dB = 10*log10(mean(nmse_sg));

  % soft decoding
  [LLR_cc_decpost_pcsi,success] = softdecode(y,z_true,mu_z_pcsi,sig2v*ones(N,1),p_s_to_y,...
        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);
  LLR_ii_decpost_pcsi = LLR_cc_decpost_pcsi([Mc-Mi+1:Mc]);% decoder info-bit estimates
  ber_pcsi = mean((LLR_ii_decpost_pcsi>0)~=ii_true); 	% average BER
  out.ber_pcsi = ber_pcsi;
  %
  [LLR_cc_decpost_bsg,success] = softdecode(y,z_hat_bsg,mu_z_bsg,sig2v*ones(N,1),p_s_to_y,...
        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);
  LLR_ii_decpost_bsg = LLR_cc_decpost_bsg([Mc-Mi+1:Mc]);% decoder info-bit estimates
  ber_bsg = mean((LLR_ii_decpost_bsg>0)~=ii_true); 	% average BER
  out.nmse_bsg = nmse_bsg; out.ber_bsg = ber_bsg;
  % 
  [LLR_cc_decpost_bgl,success] = softdecode(y,z_hat_bgl,mu_z_bgl,sig2v*ones(N,1),p_s_to_y,...
        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);
  LLR_ii_decpost_bgl = LLR_cc_decpost_bgl([Mc-Mi+1:Mc]);% decoder info-bit estimates
  ber_bgl = mean((LLR_ii_decpost_bgl>0)~=ii_true);	% average BER
  out.nmse_bgl = nmse_bgl; out.ber_bgl = ber_bgl;
  % 
  if Np>0,
    [LLR_cc_decpost_sg,success] = softdecode(y,z_hat_sg,mu_z_sg,sig2v*ones(N,1),p_s_to_y,...
        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);
    LLR_ii_decpost_sg = LLR_cc_decpost_sg([Mc-Mi+1:Mc]);% decoder info-bit estimates
    ber_sg = mean((LLR_ii_decpost_sg>0)~=ii_true);	% average BER
    out.nmse_sg = nmse_sg; out.ber_sg = ber_sg;
    % 
%    [LLR_cc_decpost_lin,success] = softdecode(y,z_hat_lin,mu_z_lin,sig2v*ones(N,1),p_s_to_y,...
%        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
%	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
%		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);
%    LLR_ii_decpost_lin = LLR_cc_decpost_lin([Mc-Mi+1:Mc]);% decoder info-bit estimates
%    ber_lin = mean((LLR_ii_decpost_lin>0)~=ii_true); 	% average BER
%    out.nmse_lin = nmse_lin; out.ber_lin = ber_lin;
  end;
 end;%genie_est


 %%%%%%%%%%%%%%%%%%
 % one-shot LASSO %
 %%%%%%%%%%%%%%%%%%
 if oneshot_lasso&&(Np>0),			% try decoding with CS-based estimate

  % bits -> symbols
  p_s_to_y = cell(1,T);			
  LLR_bb_to_del = zeros(Mb,1);			% initialize bit LLRs 
  LLR_bb_to_del(mm_train_) = -log(eps)*(2*bb_train-1); % training bits
  for t=1:T,
    p_b_to_del = 1./(1+exp(-LLR_bb_to_del((t-1)*M*Nd+[1:Nd*M]))); % LLRs -> bit probs
    p_s_to_y{t} = zeros(N,2^M);
    for i=0:Np-1,				% pilot symbol probs
      p_s_to_y{t}(n_pilot_(i+1),k_pilot_(i+1)) = 1;	 
    end;
    for i=0:Nd-1,				% data/training symbol probs
      p_Bi = ones(2^M,1)*p_b_to_del(i*M+[1:M]).';
      p_s_to_y{t}(n_data_(i+1),:) = prod(B_hyp.*p_Bi+(1-B_hyp).*(1-p_Bi),2).';
    end;
  end;%t

  % compressed channel sensing via LASSO
  x_hat_lasso = cell(1,T);
  z_hat_lasso = cell(1,T);
  mu_z_lasso = cell(1,T);
  nmse_lasso = NaN*ones(T,1);
  for t=1:T,
    A = diag(s_true{t}(n_pilot_))*Phi(n_pilot_,:);		
    nmse_lasso_test = NaN*ones(1,length(scale));
    x_hat_lasso_test = [];
    for v=1:length(scale),	% find nmse-optimal error constraint (cheating!)
      x_hat_lasso_test = spg_bpdn(A,y{t}(n_pilot_),sqrt(scale(v)*sig2v*Np),spgl1_opts);
      nmse_lasso_test(v) = norm(x_hat_lasso_test-x_true{t})^2/norm(x_true{t})^2;
    end;
    [dum,indx] = min(nmse_lasso_test);
    x_hat_lasso{t} = spg_bpdn(A,y{t}(n_pilot_),sqrt(scale(indx)*sig2v*Np),spgl1_opts);
    nmse_lasso(t) = norm(x_hat_lasso{t}-x_true{t})^2/norm(x_true{t})^2;
    z_hat_lasso{t} = Phi*x_hat_lasso{t}; 		
    mu_z_lasso{t} = sig2v*ones(N,1);		
  end;%t

  % soft decoding
  [LLR_cc_decpost,success] = softdecode(y,z_hat_lasso,mu_z_lasso,sig2v*ones(N,1),p_s_to_y,...
        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);
  LLR_ii_decpost = LLR_cc_decpost([Mc-Mi+1:Mc]);	% decoder info-bit estimates
  ber_lasso = mean((LLR_ii_decpost>0)~=ii_true);		% average BER
  out.ber_lasso = ber_lasso;
  out.nmse_lasso = nmse_lasso;
 end;%if oneshot_lasso


 %%%%%%%%%%%%%%
 % turbo GAMP %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if turbo_amp==1,
  % initialize prior to first turbo iteration 
  x_hat_amp = cell(1,T);
  z_hat_amp = cell(1,T);
  mu_z_amp = cell(1,T);
  mu_x_amp = cell(1,T);
  u_hat = cell(1,T);				% amp only
  llr_mc_e = cell(1,T);				% amp only
  llr_mc_e_old = cell(1,T);			% amp only
  llr_amp_e = cell(1,T);			% amp only
  for t=1:T,
    x_hat_amp{t} = zeros(Lx,1);		% initialize chan means
    mu_x_amp{t} = mu_onoff;			% initialize chan vars
    u_hat{t} = zeros(N,1);			% initialize normalized chan errors
    llr_mc_e{t} = -log(1./lam-1); 		% initialize MC extrinsic llrs
    llr_mc_e{t} = max(-llrmax_mc,min(llrmax_mc,llr_mc_e{t}));		% clip
  end;
  nmse_amp = NaN*ones(T,iter_turbo);
  ber_amp = NaN*ones(1,iter_turbo);
  time_amp = NaN*ones(1,iter_turbo);
  p_s_to_y = cell(1,T);
  LLR_bb_decpost = NaN*ones(Mb,1);		
  LLR_bb_decpost(mm_train_) = -log(eps)*(2*bb_train-1);% training bits
  LLR_bb_to_del = zeros(Mb,1);			% initialize bit LLRs 
  LLR_bb_to_del(mm_train_) = -log(eps)*(2*bb_train-1); % training bits
  %p_s_from_y = cell(1,T);			
  %s_hard = cell(1,T);			
  %LLR_cc_from_del = NaN*ones(Mc,1);		
  %LLR_bb_from_del = NaN*ones(Mb,1);		
  %LLR_bb_eqpost = NaN*ones(Mb,1);		

  % turbo iterations
  for turbo = 1:iter_turbo,
   tamp = tic;					% start timer
   if genie_debug, LLR_bb_to_del = -log(eps)*(2*bb_true-1); end; % perfectly known bits

   % bits -> symbols
   for t=1:T,
     p_b_to_del = 1./(1+exp(-LLR_bb_to_del((t-1)*M*Nd+[1:Nd*M]))); % LLRs -> bit probs
     p_s_to_y{t} = zeros(N,2^M);
     for i=0:Np-1,				% pilot symbol probs
       p_s_to_y{t}(n_pilot_(i+1),k_pilot_(i+1)) = 1;	 
     end;
     for i=0:Nd-1,				% data/training symbol probs
       p_Bi = ones(2^M,1)*p_b_to_del(i*M+[1:M]).';
       p_s_to_y{t}(n_data_(i+1),:) = prod(B_hyp.*p_Bi+(1-B_hyp).*(1-p_Bi),2).';
     end;
   end;%t

   % GAMP & MARKOV CHAIN DECODING
   for mc=1:max(iter_mc,1),
    for t=1:T,
     % GAMP
     for amp=1:iter_amp,
       if isnan(mu_x_amp{t}), 
         display('GAMP misconverged!'); 
	 x_hat_amp{t} = zeros(Lx,1);
	 mu_x_amp{t} = mu_onoff;
	 break; 
       end;
       % output nodes
       mu_z_amp{t} = sum(mu_x_amp{t})*ones(N,1);% mu_z_amp is constant, but still use vector 
       z_hat_amp{t} = fft(x_hat_amp{t},N);
       p_hat = z_hat_amp{t} - mu_z_amp{t}.*u_hat{t};
       zeta = p_s_to_y{t}.*exp( -(abs(y{t}*ones(1,2^M)-p_hat*s_hyp).^2 ...
  				)./(mu_z_amp{t}*s_hyp2+mu_v*ones(1,2^M))...
 		)./(mu_z_amp{t}*s_hyp2+mu_v*ones(1,2^M)) + eps; % eps prevents zeros!
       zeta = zeta./(sum(zeta,2)*ones(1,2^M));
       beta = 1./(1+1./((mu_z_amp{t}./mu_v)*s_hyp2));
       es_hat = (y{t}*s_hyp_inv-p_hat*ones(1,2^M)).*beta;
       e_hat = sum(es_hat.*zeta,2);
       mu_e = sum(( abs(es_hat-e_hat*ones(1,2^M)).^2 ...
       	+ beta.*(mu_v*abs(s_hyp_inv).^2) ).*zeta,2);
       mu_e = max(min(mu_e,ratio_ez*mu_z_amp{t}),-ratio_ez*mu_z_amp{t}); % needed modification!
       mu_u = (1-mu_e./mu_z_amp{t})./mu_z_amp{t};
       u_hat{t} = e_hat./mu_z_amp{t};
       % input nodes
       mu_q = 1/sum(mu_u);
       tmp = N*ifft(u_hat{t});
       q_hat = x_hat_amp{t} + mu_q*tmp(1:Lx);
       %nu = mu_q*(mu_on./(mu_q+mu_on));				% old method
       %gam = (q_hat.*nu)/mu_q;						% old method
       %alf = 1 + (1-lam)./lam.*mu_on./nu.*exp(-(abs(gam).^2)./nu);	% old method
       %mu_x_amp{t} = ((1-1./alf).*abs(gam).^2+nu)./alf;		% old method
       llr_amp_e{t} = -(abs(q_hat).^2).*((mu_off-mu_on ...		% support llr
       		)./(mu_on+mu_q)./(mu_off+mu_q)) - log(mu_on+mu_q) + log(mu_off+mu_q);
       llr_amp_e{t} = max(-llrmax_eq,min(llrmax_eq,llr_amp_e{t}));	% clip
       c_on = 1./(1+exp(-llr_amp_e{t}-llr_mc_e{t}) + eps);		% support posterior
       c_off = 1-c_on;
       gam_on = 1./(1+mu_q./mu_on);
       gam_off = 1./(1+mu_q./mu_off);
       %mu_x_amp{t} = c_on.*gam_on.*(mu_q+gam_on.*abs(q_hat).^2) ...
       %		+ c_off.*gam_off.*(mu_q+gam_off.*abs(q_hat).^2) ...
       %		- (abs(q_hat).^2).*(abs(c_on.*gam_on+c_off.*gam_off).^2);
       mu_x_amp{t} = c_on.*c_off.*((gam_on-gam_off).^2).*(abs(q_hat).^2) ...
       		+ (c_on.*gam_on + c_off.*gam_off).*mu_q;
       x_hat_amp_old = x_hat_amp{t};
       %x_hat_amp{t} = gam./alf;						% old method
       x_hat_amp{t} = (c_on.*gam_on+c_off.*gam_off).*q_hat;
       nmse_amp(t,turbo) = norm(x_hat_amp{t}-x_true{t})^2/norm(x_true{t})^2;
       % plot stuff
       if plot_amp_eq&&(t==1),
         figure(plot_amp_eq);
	   tit_str = ['GAMP: turbo-iter=',num2str(turbo),...
	  	  ', mc-iter=',num2str(mc),...
	  	  ', frame=',num2str(t),...
		  ', iter=',num2str(amp),...
	     	  ' (nmse=',num2str(10*log10(nmse_amp(t,turbo)),3),'dB'];
           if exist('nmse_bsg'),
	     tit_str = [tit_str,', genie=',num2str(10*log10(nmse_bsg(t)),3),'dB)'];
	   else
	     tit_str = [tit_str,')'];
           end;
         subplot(611);
           stem(supp_true{t},llr_mc_e{t}(supp_true{t}),'r.')
           hold on; stem(comp_true{t},llr_mc_e{t}(comp_true{t}),'b.'); hold off;
           axis([0,Lx,-llrmax_mc,llrmax_mc])
	   ylabel('llr-prior')
           title(tit_str);
 	 subplot(612);
	   errorbar([1:N],real(z_hat_amp{t}),3*sqrt(mu_z_amp{t}/2),3*sqrt(mu_z_amp{t}/2),'.')
	   hold on; handy = plot([1:N],real(z_true{t}),'r.'); hold off;
	   set(handy,'Markersize',7);
	   axis([0,N,-2,2]);
	   ylabel('Re(z)')
	 subplot(613);
	   errorbar([1:N],imag(z_hat_amp{t}),3*sqrt(mu_z_amp{t}/2),3*sqrt(mu_z_amp{t}/2),'.')
	   hold on; handy = plot([1:N],imag(z_true{t}),'r.'); hold off;
	   set(handy,'Markersize',7);
	   axis([0,N,-2,2])
	   ylabel('Im(z)')
	 subplot(614);
	   errorbar([1:Lx],real(x_hat_amp{t}),3*sqrt(mu_x_amp{t}/2),3*sqrt(mu_x_amp{t}/2),'.')
           if exist('nmse_bsg'),
	     hold on; handy = plot([1:Lx],real(x_hat_bsg{t}),'r.',[0,Lx],[0,0]); hold off;
           else
	     hold on; handy = plot([1:Lx],real(x_true{t}),'r.',[0,Lx],[0,0]); hold off;
           end;
	   set(handy,'Markersize',7);
	   axis([0,Lx,1.5*min(real(x_true{t})),1.5*max(real(x_true{t}))])
	   ylabel('Re(x)')
	 subplot(615);
	   errorbar([1:Lx],imag(x_hat_amp{t}),3*sqrt(mu_x_amp{t}/2),3*sqrt(mu_x_amp{t}/2),'.')
           if exist('nmse_bsg'),
	     hold on; handy = plot([1:Lx],imag(x_hat_bsg{t}),'r.',[0,Lx],[0,0]); hold off;
           else
	     hold on; handy = plot([1:Lx],imag(x_true{t}),'r.',[0,Lx],[0,0]); hold off;
           end;
	   set(handy,'Markersize',7);
	   axis([0,Lx,1.5*min(imag(x_true{t})),1.5*max(imag(x_true{t}))])
	   ylabel('Im(x)')
	 subplot(616);
	   llr_amp_p = llr_amp_e{t} + llr_mc_e{t};
           stem(supp_true{t},llr_amp_p(supp_true{t}),'r.')
           hold on; stem(comp_true{t},llr_amp_p(comp_true{t}),'b.'); hold off;
           axis([0,Lx,-llrmax_eq-llrmax_mc,llrmax_eq+llrmax_mc])
	   ylabel('llr-post')
	 drawnow;
	 %if amp==iter_amp, pause; end;
       end;%plot_amp_eq
       % break if insufficient change in channel estimate
       if norm(x_hat_amp{t}-x_hat_amp_old)^2/Lx < amp_diff_stop, break; end;
     end;%amp=1:iter_amp

     % MARKOV CHAIN DECODING
     llr_mc_e_old{t} = llr_mc_e{t};
     if iter_mc,	% MC decode if iter_mc>0
       [llr_mc_e{t},llr_mc_p] = mc_decode2(llr_amp_e{t},lam,p10,p01);	
       llr_mc_p = max(-2*llrmax_mc,min(2*llrmax_mc,llr_mc_p));		% clip
       %if t==1, [llr_amp_e{t}(1:10),llr_mc_e{t}(1:10)], end;
     end;
     llr_mc_e{t} = max(-llrmax_mc,min(llrmax_mc,llr_mc_e{t}));		% clip

    end;%t=1:T

    if mean(mean(abs([llr_mc_e{:}]-[llr_mc_e_old{:}]))) < llr_diff_stop, break; end;
   end;%mc=1:iter_mc

   %% equalizer outputs -> symbol probabilities
   %for t=1:T,
   %  % soft symbol estimates 
   %  p_s_from_y{t} = exp( -(abs(y{t}*ones(1,2^M)-z_hat_amp{t}*s_hyp).^2 ...
   %                 )./(mu_z_amp{t}*s_hyp2+mu_v*ones(1,2^M))...
   %              )./(mu_z_amp{t}*s_hyp2+mu_v*ones(1,2^M));     % ~likelihoods
   %  p_s_from_y{t} = p_s_from_y{t}./(sum(p_s_from_y{t},2)*ones(1,2^M));
   %  [dum,k_hard] = max(p_s_from_y{t}.');		% optional!
   %  s_hard{t} = s_hyp(k_hard).';			% optional!
   %end;%t
   %
   %% symbol probabilities -> bit probabilities
   %for t=1:T,
   %  p_s_tofrom_y = p_s_to_y{t}.*p_s_from_y{t};
   %  p_b_one = NaN*ones(M*Nd,1);				% proportional to Pr{b=1}
   %  p_b_zero = NaN*ones(M*Nd,1);			% proportional to Pr{b=0}
   %  for m=1:M,
   %    p_b_one([0:M:M*Nd-1]+m) = sum(p_s_tofrom_y(n_data_,k_one(m,:)),2);
   %    p_b_zero([0:M:M*Nd-1]+m) = sum(p_s_tofrom_y(n_data_,k_zero(m,:)),2);
   %  end;
   %  LLR_bb_eqpost((t-1)*M*Nd+[1:Nd*M]) = log((p_b_one+eps)./(p_b_zero+eps));
   %end;%t
   %LLR_bb_from_del = LLR_bb_eqpost - LLR_bb_to_del;
   %LLR_bb_from_del = min(max(LLR_bb_from_del,-LLRmax_eq),LLRmax_eq);
   %% I TRIED LIMITING THE LLR_bb_eqpost BUT THIS DIDNT SEEM TO WORK WELL!
   %if ~isreal(LLR_bb_from_del), error('LLR_bb_from_del is complex'); end;
   %
   %% refresh training bits
   %LLR_bb_from_del(mm_train_) = LLR_bb_to_del(mm_train_);
   %
   % LDPC decoding
   %LLR_cc_from_del(intlv) = LLR_bb_from_del(mm_data_);	% deinterleave coded bits
   %p_cc_from_del = 1./(1+exp(-LLR_cc_from_del));
   %[LLR_cc_decpost_test,success_test] = ldpc_decode1(0,1-p_cc_from_del,p_cc_from_del,HH,...
   %	LLRmax_ldpc,iter_ldpc);				% decode

   [LLR_cc_decpost,success,LLR_bb_from_del,LLR_bb_eqpost] = softdecode(...
   		y,z_hat_amp,mu_z_amp,mu_v,p_s_to_y,...
        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
   	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
   		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);
   % test to make sure "softdecode" is doing same as above
   %sum(abs(LLR_cc_decpost_test-LLR_cc_decpost))

   if ~isreal(LLR_cc_decpost),
     display(['warning: LLR_cc_decpost has ',...
     	num2str(sum(imag(LLR_cc_decpost)~=0)),' complex entries, such as from log(-eps)!'])
     %LLR_cc_decpost = real(LLR_cc_decpost);
     bad = (imag(LLR_cc_decpost)~=0);
     LLR_cc_decpost(bad) = zeros(sum(bad),1);
   end;
   LLR_ii_decpost = LLR_cc_decpost([Mc-Mi+1:Mc]);	% decoder info-bit estimates
   LLR_bb_decpost(mm_data_) = LLR_cc_decpost(intlv);	% reinterleave
   LLR_bb_to_del_old = LLR_bb_to_del;			% previous decoder extrinsic outs 
   % decoder extrinsic outs
   LLR_bb_to_del(mm_data_) = LLR_bb_decpost(mm_data_)-LLR_bb_from_del(mm_data_);	
   LLR_diff = mean(abs(LLR_bb_to_del-LLR_bb_to_del_old));  % mean LLR difference

   % examine output
   ber_amp(turbo) = mean((LLR_ii_decpost>0)~=ii_true);	% average info-bit BER
   %ber_amp(turbo) = mean((LLR_bb_decpost>0)~=bb_true);	% average coded-bit BER

   % plot LLRs
   if plot_amp_dec,
     figure(plot_amp_dec);
     tit_str = ['GAMP: turbo iteration ',num2str(turbo),'  (ber=',num2str(ber_amp(turbo),5),...
      	', nmse=',num2str(10*log10(mean(nmse_amp(:,turbo))),3),'dB'];
     if exist('ber_pcsi')
       tit_str = [tit_str,', ber-pcsi=',num2str(ber_pcsi,5),')'];
     else
       tit_str = [tit_str,')'];
     end;
     subplot(411)
       plot(indx_bb1,LLR_bb_to_del_old(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_to_del_old(indx_bb0),'b.',...
 	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('eq in')
       title(tit_str);
     subplot(412)
       plot(indx_bb1,LLR_bb_eqpost(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_eqpost(indx_bb0),'b.',...
 	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('eq post')
     subplot(413)
       plot(indx_bb1,LLR_bb_decpost(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_decpost(indx_bb0),'b.',...
	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('dec post')
     subplot(414)
       plot(indx_bb1,LLR_bb_to_del(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_to_del(indx_bb0),'b.',...
	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('dec out')
     drawnow;
   end;
   %if turbo~=iter_turbo, pause; end;

   time_amp(turbo) = toc(tamp);
   if stop_ber==1,
     if ber_amp(turbo)==0, break; end;		% stop after zero BER 
   elseif stop_ber==2,
     if success, break; end;			% stop after ldpc decoder indicates success
   end;
   if LLR_diff<LLR_diff_stop, break; end;
   if genie_debug; break; end;
  end;%turbo

  % return outputs
  out.ber_amp = ber_amp;
  out.nmse_amp = nmse_amp;
  out.time_amp = time_amp;
 end;%turbo_amp


 %%%%%%%%%%%%%%%
 % turbo LMMSE %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if turbo_lin==1,
  % initialize prior to first turbo iteration 
  x_hat_lin = cell(1,T);
  z_hat_lin = cell(1,T);
  mu_z_lin = cell(1,T);
  mu_x_lin = cell(1,T);	
  nmse_lin = NaN*ones(T,iter_turbo);
  ber_lin = NaN*ones(1,iter_turbo);
  time_lin = NaN*ones(1,iter_turbo);
  p_s_to_y = cell(1,T);
  LLR_bb_decpost = NaN*ones(Mb,1);		
  LLR_bb_decpost(mm_train_) = -log(eps)*(2*bb_train-1);% training bits
  LLR_bb_to_del = zeros(Mb,1);			% initialize bit LLRs 
  LLR_bb_to_del(mm_train_) = -log(eps)*(2*bb_train-1); % training bits

  % turbo iterations
  for turbo = 1:iter_turbo,
   tlin = tic;					% start timer
   if genie_debug, LLR_bb_to_del = -log(eps)*(2*bb_true-1); end; % perfectly known bits

   % bits -> symbols
   for t=1:T,
     p_b_to_del = 1./(1+exp(-LLR_bb_to_del((t-1)*M*Nd+[1:Nd*M]))); % LLRs -> bit probs
     p_s_to_y{t} = zeros(N,2^M);
     for i=0:Np-1,				% pilot symbol probs
       p_s_to_y{t}(n_pilot_(i+1),k_pilot_(i+1)) = 1;	 
     end;
     for i=0:Nd-1,				% data/training symbol probs
       p_Bi = ones(2^M,1)*p_b_to_del(i*M+[1:M]).';
       p_s_to_y{t}(n_data_(i+1),:) = prod(B_hyp.*p_Bi+(1-B_hyp).*(1-p_Bi),2).';
     end;
   end;%t

   % SOFT LMMSE CHANNEL ESTIMATION
   for t=1:T,
     s_hat = p_s_to_y{t}*s_hyp.';			% soft symbol mean
     mu_s = p_s_to_y{t}*s_hyp2.' - abs(s_hat).^2;	% soft symbol variances
     if hard_fdbk,
       [dum,k_hard] = max(p_s_to_y{t}.');
       s_hat = s_hyp(k_hard).';				% hard symbol est
     end;
     
     A = diag(s_hat)*Phi; 
     Rxy = diag(mu_onoff)*A';
     %Ryy = A*diag(mu_onoff)*A'+diag(sig2v+mu_s.*diag(Phi*diag(mu_onoff)*Phi'));
     Ryy = A*diag(mu_onoff)*A'+diag(sig2v+mu_s*sum(mu_onoff));%shortcut
     %RxyRyy_inv = Rxy*inv(Ryy);
     RxyRyy_inv = Rxy/Ryy;
     x_hat_lin{t} = RxyRyy_inv*y{t};			% LMMSE chan estimate
     %Mu_x_lin = diag(mu_onoff) - RxyRyy_inv*Rxy';	% LMMSE covariance matrix
     %mu_x_lin{t} = abs(diag(Mu_x_lin));		% LMMSE variances
     mu_x_lin{t} = mu_onoff - abs(sum(RxyRyy_inv.*conj(Rxy),2));% LMMSE variances shortcut
     z_hat_lin{t} = fft(x_hat_lin{t},N);
     %mu_z_lin{t} = abs(diag(Phi*Mu_x_lin*Phi'));	% exact covariance
     %mu_z_lin{t} = abs(diag(Phi*diag(mu_x_lin{t})*Phi'));% approximate covariance
     mu_z_lin{t} = sum(mu_x_lin{t})*ones(N,1);% approximate covariance shortcut
     nmse_lin(t,turbo) = norm(x_hat_lin{t}-x_true{t})^2/norm(x_true{t})^2;

     % plot stuff
     if plot_lin_eq&&(t==1),
       figure(plot_lin_eq);
	   tit_str = ['LMMSE: turbo-iter=',num2str(turbo),...
	  	  ', frame=',num2str(t),...
	     	  ' (nmse=',num2str(10*log10(nmse_lin(t,turbo)),3),'dB'];
           if exist('nmse_bsg'),
	     tit_str = [tit_str,', genie=',num2str(10*log10(nmse_bsg(t)),3),'dB)'];
	   else
	     tit_str = [tit_str,')'];
           end;
 	 subplot(411);
	   errorbar([1:N],real(z_hat_lin{t}),3*sqrt(mu_z_lin{t}/2),3*sqrt(mu_z_lin{t}/2),'.')
	   hold on; handy = plot([1:N],real(z_true{t}),'r.'); hold off;
	   set(handy,'Markersize',7);
	   axis([0,N,-2,2]);
	   ylabel('Re(z)')
           title(tit_str);
	 subplot(412);
	   errorbar([1:N],imag(z_hat_lin{t}),3*sqrt(mu_z_lin{t}/2),3*sqrt(mu_z_lin{t}/2),'.')
	   hold on; handy = plot([1:N],imag(z_true{t}),'r.'); hold off;
	   set(handy,'Markersize',7);
	   axis([0,N,-2,2])
	   ylabel('Im(z)')
	 subplot(413);
	   errorbar([1:Lx],real(x_hat_lin{t}),3*sqrt(mu_x_lin{t}/2),3*sqrt(mu_x_lin{t}/2),'.')
           if exist('nmse_bsg'),
	     hold on; handy = plot([1:Lx],real(x_hat_bsg{t}),'r.',[0,Lx],[0,0]); hold off;
           else
	     hold on; handy = plot([1:Lx],real(x_true{t}),'r.',[0,Lx],[0,0]); hold off;
           end;
	   set(handy,'Markersize',7);
	   axis([0,Lx,1.5*min(real(x_true{t})),1.5*max(real(x_true{t}))])
	   ylabel('Re(x)')
	 subplot(414);
	   errorbar([1:Lx],imag(x_hat_lin{t}),3*sqrt(mu_x_lin{t}/2),3*sqrt(mu_x_lin{t}/2),'.')
           if exist('nmse_bsg'),
	     hold on; handy = plot([1:Lx],imag(x_hat_bsg{t}),'r.',[0,Lx],[0,0]); hold off;
           else
	     hold on; handy = plot([1:Lx],imag(x_true{t}),'r.',[0,Lx],[0,0]); hold off;
           end;
	   set(handy,'Markersize',7);
	   axis([0,Lx,1.5*min(imag(x_true{t})),1.5*max(imag(x_true{t}))])
	   ylabel('Im(x)')
	 drawnow;
     end;%plot_lin_eq
   end;%t=1:T

   [LLR_cc_decpost,success,LLR_bb_from_del,LLR_bb_eqpost] = softdecode(...
   		y,z_hat_lin,mu_z_lin,mu_v,p_s_to_y,...
        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
   	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
   		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);

   if ~isreal(LLR_cc_decpost),
     display(['warning: LLR_cc_decpost has ',...
     	num2str(sum(imag(LLR_cc_decpost)~=0)),' complex entries, such as from log(-eps)!'])
     %LLR_cc_decpost = real(LLR_cc_decpost);
     bad = (imag(LLR_cc_decpost)~=0);
     LLR_cc_decpost(bad) = zeros(sum(bad),1);
   end;
   LLR_ii_decpost = LLR_cc_decpost([Mc-Mi+1:Mc]);	% decoder info-bit estimates
   LLR_bb_decpost(mm_data_) = LLR_cc_decpost(intlv);	% reinterleave
   LLR_bb_to_del_old = LLR_bb_to_del;			% previous decoder extrinsic outs 
   % decoder extrinsic outs
   LLR_bb_to_del(mm_data_) = LLR_bb_decpost(mm_data_)-LLR_bb_from_del(mm_data_);	
   LLR_diff = mean(abs(LLR_bb_to_del-LLR_bb_to_del_old));  % mean LLR difference

   % examine output
   ber_lin(turbo) = mean((LLR_ii_decpost>0)~=ii_true);	% average info-bit BER
   %ber_lin(turbo) = mean((LLR_bb_decpost>0)~=bb_true);	% average coded-bit BER

   % plot LLRs
   if plot_lin_dec,
     figure(plot_lin_dec);
     tit_str = ['LMMSE: turbo iteration ',num2str(turbo),'  (ber=',num2str(ber_lin(turbo),5),...
      	', nmse=',num2str(10*log10(mean(nmse_lin(:,turbo))),3),'dB'];
     if exist('ber_pcsi')
       tit_str = [tit_str,', ber-pcsi=',num2str(ber_pcsi,5),')'];
     else
       tit_str = [tit_str,')'];
     end;
     subplot(411)
       plot(indx_bb1,LLR_bb_to_del_old(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_to_del_old(indx_bb0),'b.',...
 	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('eq in')
       title(tit_str);
     subplot(412)
       plot(indx_bb1,LLR_bb_eqpost(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_eqpost(indx_bb0),'b.',...
 	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('eq post')
     subplot(413)
       plot(indx_bb1,LLR_bb_decpost(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_decpost(indx_bb0),'b.',...
	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('dec post')
     subplot(414)
       plot(indx_bb1,LLR_bb_to_del(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_to_del(indx_bb0),'b.',...
	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('dec out')
     drawnow;
   end;
   %if turbo~=iter_turbo, pause; end;

   time_lin(turbo) = toc(tlin);
   if stop_ber==1,
     if ber_lin(turbo)==0, break; end;		% stop after zero BER 
   elseif stop_ber==2,
     if success, break; end;			% stop after ldpc decoder indicates success
   end;
   if LLR_diff<LLR_diff_stop, break; end;
   if genie_debug; break; end;
  end;%turbo

  % return outputs
  out.ber_lin = ber_lin;
  out.nmse_lin = nmse_lin;
  out.time_lin = time_lin;
 end;%turbo_lin


 %%%%%%%%%%%%%%%
 % turbo LASSO %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if turbo_las==1,
  % initialize prior to first turbo iteration 
  x_hat_las = cell(1,T);
  z_hat_las = cell(1,T);
  mu_z_las = cell(1,T);
  mu_x_las = cell(1,T);	
  nmse_las = NaN*ones(T,iter_turbo);
  ber_las = NaN*ones(1,iter_turbo);
  time_las = NaN*ones(1,iter_turbo);
  p_s_to_y = cell(1,T);
  LLR_bb_decpost = NaN*ones(Mb,1);		
  LLR_bb_decpost(mm_train_) = -log(eps)*(2*bb_train-1);% training bits
  LLR_bb_to_del = zeros(Mb,1);			% initialize bit LLRs 
  LLR_bb_to_del(mm_train_) = -log(eps)*(2*bb_train-1); % training bits
  mse_las_test_min_old = inf*ones(1,iter_turbo);	% lasso MSEs

  % turbo iterations
  for turbo = 1:iter_turbo,
   tlasso = tic;				% start timer
   if genie_debug, LLR_bb_to_del = -log(eps)*(2*bb_true-1); end; % perfectly known bits

   % bits -> symbols
   for t=1:T,
     p_b_to_del = 1./(1+exp(-LLR_bb_to_del((t-1)*M*Nd+[1:Nd*M]))); % LLRs -> bit probs
     p_s_to_y{t} = zeros(N,2^M);
     for i=0:Np-1,				% pilot symbol probs
       p_s_to_y{t}(n_pilot_(i+1),k_pilot_(i+1)) = 1;	 
     end;
     for i=0:Nd-1,				% data/training symbol probs
       p_Bi = ones(2^M,1)*p_b_to_del(i*M+[1:M]).';
       p_s_to_y{t}(n_data_(i+1),:) = prod(B_hyp.*p_Bi+(1-B_hyp).*(1-p_Bi),2).';
     end;
   end;%t

   % SOFT LASSO CHANNEL ESTIMATION
   spgl1_stat = zeros(1,T);
   for t=1:T,
     if (turbo==1)&&(Np>0),
       % use only pilot observations
       A = diag(s_true{t}(n_pilot_))*Phi(n_pilot_,:);		
       mse_las_test = NaN*ones(1,length(spgl1_scale));
       X_hat_las_test = NaN*ones(Lx,length(spgl1_scale));
       for v=1:length(spgl1_scale),	% find nmse-optimal error constraint (cheating!)
         X_hat_las_test(:,v) = spg_bpdn(...
	 	A,y{t}(n_pilot_),sqrt(spgl1_scale(v)*sig2v*Np),spgl1_opts);
         mse_las_test(v) = norm(X_hat_las_test(:,v)-x_true{t})^2/Lx;
       end;
       [mse_las_test_min,indx] = min(mse_las_test);
       x_hat_las{t} = X_hat_las_test(:,indx);
       nmse_las(t,turbo) = mse_las_test_min*Lx/norm(x_true{t})^2;
       mu_x_las{t} = mse_las_test_min*ones(Lx,1);
       z_hat_las{t} = fft(x_hat_las{t},N);
       mu_z_las{t} = Lx*mse_las_test_min*ones(N,1);
       mse_las_test_min_old(t) = mse_las_test_min;
     else%if (turbo~=1)||(Np==0)
       s_hat = p_s_to_y{t}*s_hyp.';			% soft symbol mean
       mu_s = p_s_to_y{t}*s_hyp2.' - abs(s_hat).^2;	% soft symbol variances
       if hard_fdbk,
         [dum,k_hard] = max(p_s_to_y{t}.');
         s_hat = s_hyp(k_hard).';			% hard symbol est
       end;
     
       mu_n = sig2v+mu_s*sum(mu_onoff); 
       A = diag(s_hat./sqrt(mu_n))*Phi;	% whitened-channel mtx
       mse_las_test = NaN*ones(1,length(spgl1_scale));
       X_hat_las_test = NaN*ones(Lx,length(spgl1_scale));
       spgl1_stat_test = NaN*ones(1,length(spgl1_scale));
       for v=1:length(spgl1_scale),	% find nmse-optimal error constraint (cheating!)
 	 if 1,%if v==1,			% cold start
           [X_hat_las_test(:,v),resid,grad,info] = spg_bpdn(...
	   	A,y{t}./sqrt(mu_n),sqrt(spgl1_scale(v)*N),spgl1_opts);
	 else				% warm start
           [X_hat_las_test(:,v),resid,grad,info] = spgl1(...
	   	A,y{t}./sqrt(mu_n),0,sqrt(spgl1_scale(v)*N),X_hat_las_test(:,v-1),spgl1_opts);
	 end;
         mse_las_test(v) = norm(X_hat_las_test(:,v)-x_true{t})^2/Lx;
       spgl1_stat_test(v) = info.stat;
       end;
       [mse_las_test_min,indx] = min(mse_las_test);
       if mse_las_test_min < mse_las_test_min_old(t),	% if estimate has improved...
         spgl1_stat(t) = spgl1_stat_test(indx);
         x_hat_las{t} = X_hat_las_test(:,indx);
         nmse_las(t,turbo) = mse_las_test_min*Lx/norm(x_true{t})^2;
         mu_x_las{t} = mse_las_test_min*ones(Lx,1);
         z_hat_las{t} = fft(x_hat_las{t},N);
         mu_z_las{t} = Lx*mse_las_test_min*ones(N,1);
         mse_las_test_min_old(t) = mse_las_test_min;
       else,
         nmse_las(t,turbo) = nmse_las(t,turbo-1); 
       end;
     end;%if (turbo==1)&&(Np>0)

     % plot stuff
     if plot_las_eq&&(t==1),
       figure(plot_las_eq);
	   tit_str = ['LASSO: turbo-iter=',num2str(turbo),...
	  	  ', frame=',num2str(t),...
	     	  ' (nmse=',num2str(10*log10(nmse_las(t,turbo)),3),'dB'];
           if exist('nmse_bsg'),
	     tit_str = [tit_str,', genie=',num2str(10*log10(nmse_bsg(t)),3),'dB)'];
	   else
	     tit_str = [tit_str,')'];
           end;
 	 subplot(411);
	   errorbar([1:N],real(z_hat_las{t}),3*sqrt(mu_z_las{t}/2),3*sqrt(mu_z_las{t}/2),'.')
	   hold on; handy = plot([1:N],real(z_true{t}),'r.'); hold off;
	   set(handy,'Markersize',7);
	   axis([0,N,-2,2]);
	   ylabel('Re(z)')
           title(tit_str);
	 subplot(412);
	   errorbar([1:N],imag(z_hat_las{t}),3*sqrt(mu_z_las{t}/2),3*sqrt(mu_z_las{t}/2),'.')
	   hold on; handy = plot([1:N],imag(z_true{t}),'r.'); hold off;
	   set(handy,'Markersize',7);
	   axis([0,N,-2,2])
	   ylabel('Im(z)')
	 subplot(413);
	   errorbar([1:Lx],real(x_hat_las{t}),3*sqrt(mu_x_las{t}/2),3*sqrt(mu_x_las{t}/2),'.')
           if exist('nmse_bsg'),
	     hold on; handy = plot([1:Lx],real(x_hat_bsg{t}),'r.',[0,Lx],[0,0]); hold off;
           else
	     hold on; handy = plot([1:Lx],real(x_true{t}),'r.',[0,Lx],[0,0]); hold off;
           end;
	   set(handy,'Markersize',7);
	   axis([0,Lx,1.5*min(real(x_true{t})),1.5*max(real(x_true{t}))])
	   ylabel('Re(x)')
	 subplot(414);
	   errorbar([1:Lx],imag(x_hat_las{t}),3*sqrt(mu_x_las{t}/2),3*sqrt(mu_x_las{t}/2),'.')
           if exist('nmse_bsg'),
	     hold on; handy = plot([1:Lx],imag(x_hat_bsg{t}),'r.',[0,Lx],[0,0]); hold off;
           else
	     hold on; handy = plot([1:Lx],imag(x_true{t}),'r.',[0,Lx],[0,0]); hold off;
           end;
	   set(handy,'Markersize',7);
	   axis([0,Lx,1.5*min(imag(x_true{t})),1.5*max(imag(x_true{t}))])
	   ylabel('Im(x)')
	 drawnow;
     end;%plot_las_eq
   end;%t=1:T

   [LLR_cc_decpost,success,LLR_bb_from_del,LLR_bb_eqpost] = softdecode(...
   		y,z_hat_las,mu_z_las,mu_v,p_s_to_y,...
        	LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
   	        intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
   		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_);

   if ~isreal(LLR_cc_decpost),
     display(['warning: LLR_cc_decpost has ',...
     	num2str(sum(imag(LLR_cc_decpost)~=0)),' complex entries, such as from log(-eps)!'])
     %LLR_cc_decpost = real(LLR_cc_decpost);
     bad = (imag(LLR_cc_decpost)~=0);
     LLR_cc_decpost(bad) = zeros(sum(bad),1);
   end;
   LLR_ii_decpost = LLR_cc_decpost([Mc-Mi+1:Mc]);	% decoder info-bit estimates
   LLR_bb_decpost(mm_data_) = LLR_cc_decpost(intlv);	% reinterleave
   LLR_bb_to_del_old = LLR_bb_to_del;			% previous decoder extrinsic outs 
   % decoder extrinsic outs
   LLR_bb_to_del(mm_data_) = LLR_bb_decpost(mm_data_)-LLR_bb_from_del(mm_data_);	
   LLR_diff = mean(abs(LLR_bb_to_del-LLR_bb_to_del_old));  % mean LLR difference
   % examine output
   ber_las(turbo) = mean((LLR_ii_decpost>0)~=ii_true);	% average info-bit BER
   %ber_las(turbo) = mean((LLR_bb_decpost>0)~=bb_true);	% average coded-bit BER

   % plot LLRs
   if plot_las_dec,
     figure(plot_las_dec);
     tit_str = ['LASSO: turbo iteration ',num2str(turbo),'  (ber=',num2str(ber_las(turbo),5),...
      	', nmse=',num2str(10*log10(mean(nmse_las(:,turbo))),3),'dB'];
     if exist('ber_pcsi')
       tit_str = [tit_str,', ber-pcsi=',num2str(ber_pcsi,5),')'];
     else
       tit_str = [tit_str,')'];
     end;
     subplot(411)
       plot(indx_bb1,LLR_bb_to_del_old(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_to_del_old(indx_bb0),'b.',...
 	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('eq in')
       title(tit_str);
     subplot(412)
       plot(indx_bb1,LLR_bb_eqpost(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_eqpost(indx_bb0),'b.',...
 	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('eq post')
     subplot(413)
       plot(indx_bb1,LLR_bb_decpost(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_decpost(indx_bb0),'b.',...
	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('dec post')
     subplot(414)
       plot(indx_bb1,LLR_bb_to_del(indx_bb1),'r.',...
 	indx_bb0,LLR_bb_to_del(indx_bb0),'b.',...
	[0,Mc],[0,0],'k');
       axis([0,Mc,-LLRmax_ldpc-5,LLRmax_ldpc+5])
       ylabel('dec out')
     drawnow;
   end;
   %if turbo~=iter_turbo, pause; end;

   time_las(turbo) = toc(tlasso);
   if stop_ber==1,
     if ber_las(turbo)==0, break; end;		% stop after zero BER 
   elseif stop_ber==2,
     if success, break; end;			% stop after ldpc decoder indicates success
   end;
   if LLR_diff<LLR_diff_stop, break; end;
   if sum(spgl1_stat>4), display('breaking due to SPGL1 error'); break; end;
   if genie_debug; break; end;
  end;%turbo

  % return outputs
  out.ber_las = ber_las;
  out.nmse_las = nmse_las;
  out.time_las = time_las;
 end;%turbo_las
