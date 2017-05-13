 % simulation parameters
 M = 100; 		% number of channel realizations to generate [dflt=50000]
 thresh_dB = -150;	% precision of propagation model: at least -150dB!! [dflt=-200]
 plot_chan_lin = 1;	% plot amplitude realizations
 plot_chan_dB = 2;	% plot amplitude realizations
 plot_taphist = 3;	% plot histograms of selected taps
 plot_stats = 4;	% plot per-tap statistics 
 plot_em = 5;		% plot per-tap EM stuff 
 plot_mc = 6;		% plot per-tap markov-chain stuff 
 plot_states = 7;	% plot state realizations 
 plot_corr = 8;		% plot correlations
 print_plots = 0;	% prints plots to .eps
 save_stats = 0; 

 % OFDM parameters 
 L = 256;		% total number of channel taps [dflt=256]
 T = 1.0e-6/L;		% baud interval [dflt=1.0e-6/L]
 Lpre = 20;		% number of channel taps before first arrival [dflt=20]
 Lam0 = 1/T;		% first cluster arrival rate (Hz)...[1/T or 2/T]
 alf_rc = 0.5;		% raised-cosine parameter [dflt=0.5]

 % set propagation channel parameters
 prop.chan_type = 2;	% channel type: 1=outdoor-LOS, 2=outdoor-NLOS
 prop.thresh_dB = -150;	% precision of propagation model: [dflt=-150]
 prop.Lp = 100;		% number of paths per cluster [dflt=100]
 prop.Lam0 = Lam0;	% first cluster arrival rate (Hz)

 % modeling parameters
 iid_bg = 0;		% restrict model to iid bernoulli-gaussian [dflt=0]

 % set random seeds
 % set randn and rand seeds
 reset(RandStream.getDefaultStream);   % uncomment to start at reference
 %defaultStream.State = savedState;     % uncomment to use previous state
 defaultStream = RandStream.getDefaultStream;   % leave uncommented!
 savedState = defaultStream.State;      % leave uncommented!


% test with gaussian mixture 
%test = 0;
%if test,
% h = zeros(M,L);
% lam_true = 0.39;
% var_on_true = 0.0064;
% var_off_true = 4e-14;
% K = round(L*lam_true);
% h(:,1:K) = sqrt(var_off_true)*(randn(M,K)+j*randn(M,K))/sqrt(2);
% h(:,K+1:L) = sqrt(var_on_true)*(randn(M,L-K)+j*randn(M,L-K))/sqrt(2);
% E_h = (norm(h,'fro')^2)/M; 
% var_on_true = var_on_true/E_h 
% var_off_true = var_off_true/E_h 
%end;


if 1,
 % generate propagation channels and ofdm channels
 n = [0:L-1]+eps;		% baud-normalized sampling times
 h = zeros(M,L);		% ofdm channels
 fprintf(1,'percent complete:  0')
 for i=1:M,
   fprintf(1,'\b\b%2d',floor((i-1)/M*100));
   out = prop_chan(prop);	% generate propagation {delays,amplitudes}
   %out.tau = i/M*T; out.amp = 1;	% trivial channels for testing
   for p=1:length(out.tau),	% for each path... 
     n_p = Lpre+out.tau(p)/T;	% baud-normalized propagation delay
     h(i,:) = h(i,:) + exp(j*2*pi*rand(1))*out.amp(p)*cos(...
	pi*alf_rc*(n-n_p))./(1-(2*alf_rc*(n-n_p)).^2).*sinc(n-n_p);
   end;
   tau_max = max(out.tau)+Lpre*T;

   num_plot = 6;
   if plot_chan_dB&(i<num_plot+1),
     figure(plot_chan_dB);
     subplot(2,num_plot/2,i)
       handy = stem(Lpre+out.tau/T,10*log10(out.amp),'r'); 
       set(handy,'MarkerSize',0);
       hold on;
         plot([1:L],10*log10(abs(h(i,:))),'.-');
       hold off;
       axis([0,max(tau_max/T,L),thresh_dB+50,5]);
       grid on;
       drawnow;
   end;%plot_chan_lin
   if plot_chan_lin&(i<num_plot+1),
     figure(plot_chan_lin); 
     subplot(2,num_plot/2,i)
       handy = stem(Lpre+out.tau/T,out.amp,'r'); set(handy,'MarkerSize',0);
       hold on;
         plot([1:L],abs(h(i,:)),'.-');
       hold off;
       axe = axis; axis([0,max(tau_max/T,L),0,axe(4)]);
       grid on;
       drawnow;
   end;%plot_chan_dB
 end;%M
 fprintf(1,'\b\b100\n')
end

 % normalize and analyze ofdm channels
 E_h = (norm(h,'fro')^2)/M; 		% compute average channel energy
 h = h/sqrt(E_h);			% normalize channels!
 mu_h = mean(h);			% per-tap mean
 var_h = var(h);			% per-tap variance
 pdp_h = abs(mu_h).^2 + var_h;		% power delay profile
 cm4_h = mean(abs(h-ones(M,1)*mu_h).^4);% 4th central moment
 kurt_h = cm4_h./(var_h.^2);		% kurtosis
 ell1_h = sum(abs(h),2);		% ell-1 norm
 ell1_h_norm = sum(abs(diag(1./(sqrt(diag(h*h'))))*h),2); % ell-1 norm after ell-2 calib
 if plot_stats,
   figure(plot_stats); 
   subplot(211)
     plot([1:L],10*log10(pdp_h));
     grid on;
     title('power delay profile')
     axe=axis; axis([0,L,thresh_dB,axe(4)]);
   subplot(212)
     semilogy([1:L],kurt_h)
     grid on;
     title('kurtosis')
     axe=axis; axis([0,L,axe(3:4)]);
     xlabel('lag [baud]');
 end;%stats

 % histograms
 if plot_taphist,
   figure(plot_taphist)
   bins = 50;
   subplot(221);
     tap=round(Lpre/4);
%tap=round(Lpre+3);
     hist(real(h(:,tap)),bins); 
     %hist(20*log10(abs(h(:,tap))),bins); 
     grid on;
     axis('tight')
     title(['histogram at lag ',num2str(tap)])
   subplot(222);
     tap=round(Lpre+3);
%tap=round(0.25*L);
     hist(real(h(:,tap)),bins); 
     %hist(20*log10(abs(h(:,tap))),bins); 
     grid on;
     axis('tight')
     title(['histogram at lag ',num2str(tap)])
   subplot(223);
     tap=round(0.5*L);
%tap=round(0.5*L);
     hist(real(h(:,tap)),bins); 
     %hist(20*log10(abs(h(:,tap))),bins); 
     grid on;
     axis('tight')
     title(['histogram at lag ',num2str(tap)])
   subplot(224);
     tap=round(0.9*L);
%tap=round(0.75*L);
     hist(real(h(:,tap)),bins); 
     %hist(20*log10(abs(h(:,tap))),bins); 
     grid on;
     axis('tight')
     title(['histogram at lag ',num2str(tap)])
 end;%taphist

 % fit 2-state gaussian-mixture model using EM
 iter_max = 50;
 diff_lam = 1e-4;			% [dflt 1e-4]
 mu_gm = zeros(2,L);			% initial GM means {off;on}
 if iid_bg,
   var_gm = mean(pdp_h)*[0.1;2]*ones(1,L);
   lam_gm = 0.5*ones(1,L);
 else
   var_gm = [0.4*pdp_h(1:Lpre),0.4*pdp_h(Lpre+1:L);...	% initial GM variances
 	   0.6*pdp_h(1:Lpre),0.6*pdp_h(Lpre+1:L)];	
   lam_gm = [0.1*ones(1,1.5*Lpre),0.7*ones(1,L-1.5*Lpre)];% initial GM off-probs
 end;
 for i=1:iter_max,
   Voff = ones(M,1)*var_gm(1,:);
   Von = ones(M,1)*var_gm(2,:);
   CNoff = exp(-(abs(h-ones(M,1)*mu_gm(1,:)).^2)./Voff)./Voff;
   CNon = exp(-(abs(h-ones(M,1)*mu_gm(2,:)).^2)./Von)./Von;
   Lamoff = ones(M,1)*lam_gm;
   Lamon = 1-Lamoff;
   gamoff = (Lamoff.*CNoff)./(Lamoff.*CNoff+Lamon.*CNon+eps);
   gamon = 1-gamoff;				% support posterior
   %mu_gm(1,:) = sum(gamoff.*h,1)./Moff;	% comment out for zero-mean
   %mu_gm(2,:) = sum(gamon.*h,1)./Mon;		% comment out for zero-mean
   if iid_bg,					% iid Bernoulli-Gaussian
     MLoff = sum(gamoff(:),1);
     MLon = M*L-MLoff;
     var_gm(1,:) = sum(sum(gamoff.*(abs(h-ones(M,1)*mu_gm(1,:)).^2)))/MLoff*ones(1,L);
     var_gm(2,:) = sum(sum(gamon.*(abs(h-ones(M,1)*mu_gm(2,:)).^2)))/MLon*ones(1,L);
     lam_gm_old = lam_gm;
     lam_gm = MLoff/M/L*ones(1,L);
   else						% nid 2-state Gaussian mixture
     Moff = sum(gamoff,1);
     Mon = M-Moff;
     var_gm(1,:) = sum(gamoff.*(abs(h-ones(M,1)*mu_gm(1,:)).^2))./Moff;
     var_gm(2,:) = sum(gamon.*(abs(h-ones(M,1)*mu_gm(2,:)).^2))./Mon;
     lam_gm_old = lam_gm;
     lam_gm = Moff/M;
   end;
%figure(8);
%imagesc(gamon); colorbar;
%drawnow
%[lam_gm(1), var_gm(:,1)']
   if mean(abs(lam_gm-lam_gm_old))<diff_lam, break; end;

   if plot_em,
     figure(plot_em); 
     subplot(311)
       plot([1:L],1-lam_gm);
       grid on;
       ylabel('lambda')
       axis([0,L,0,1]);
       title(['gaussian-mixture fitting, EM iter=',num2str(i)])
     subplot(312)
       handy = plot([1:L],10*log10(var_gm(2,:)./pdp_h),'r',...
	 	[1:L],10*log10(var_gm(1,:)./pdp_h));
       grid on;
       ylabel('variance/PDP [dB]')
       axis('tight');%axe=axis; axis([0,L,axe(3:4)]);
       legend(handy,'on','off','Location','SouthEast')
     %subplot(313)
       %handy = plot([1:L],abs(mu_gm(1,:))./sqrt(pdp_h),'r',...
       %	 	[1:L],abs(mu_gm(2,:))./sqrt(pdp_h),'--');
       %ylabel('|mean|/sqrt(dpp)')
       %axe=axis; axis([1,L,axe(3:4)]);
       %legend(handy,'off','on')
     subplot(313)
       plot([1:L],10*log10(pdp_h));
       grid on;
       ylabel('PDP [dB]')
       axis('tight');%axe=axis; axis([0,L,thresh_dB,axe(4)]);
       xlabel('lag [baud]');
     drawnow;
   end;
 end;%i

 % fit Markov model to sparsity pattern
 tapon = double(gamon>0.5);
 tapoff = 1-tapon;
 if iid_bg,
   p1 = mean(tapon(:))*ones(1,L); p0 = 1-p1;
   p10 = p1(1:L-1);
   p01 = p0(1:L-1);
 else
   p1 = mean(tapon); p0 = 1-p1;
   p01 = mean(tapon(:,1:L-1).*tapoff(:,2:L)+eps)./(p1(1:L-1)+eps); 
   p10 = mean(tapoff(:,1:L-1).*tapon(:,2:L)+eps)./(p0(1:L-1)+eps); 
 end;

 % save channel statistics
 if save_stats,
   if iid_bg, 
     bggm_str = 'bgstats_type='; 
   else 
     bggm_str = 'gmstats_type='; 
   end;
   save_str = [bggm_str,num2str(prop.chan_type),'_L=',num2str(L),'_TL=',num2str(T*L*1e6),'_alf=',num2str(alf_rc),'_TLam=',num2str(T*Lam0),'.mat'];
   mu_on = var_gm(2,:)';
   mu_off = var_gm(1,:)';
   lam = p1';
   save(save_str,'lam','p01','p10','mu_on','mu_off','E_h');
 end;

 % simulate markov chain
 Msim = min(5000,M);			% number of realizations to simulate
 son = NaN*ones(Msim,L); 
 for m=1:Msim,
   son(m,1) = ( rand(1)<p1(1) );	% initial state
   for l=1:L-1, 			% markov chain realization
     son(m,l+1) = (rand(1) < [son(m,l),1-son(m,l)]*[1-p01(l);p10(l)]); 
   end;
 end;
 soff = 1-son;
 p1_test = mean(son); p0_test = 1-p1_test;
 p01_test = mean(son(:,1:L-1).*soff(:,2:L)+eps)./(p1_test(:,1:L-1)+eps);
 p10_test = mean(soff(:,1:L-1).*son(:,2:L)+eps)./(p0_test(:,1:L-1)+eps);
 if plot_mc,
  figure(plot_mc); 
  subplot(211)
    handy = plot([1:L-1],p01,'r',[1:L-1],p01_test);
    %handy = plot([1:L-1],1./(1+p01./p10),'r',[1:L-1],1./(1+p01_test./p10_test));
    axe=axis; axis([1,L,axe(3:4)]);
    legend(handy,'measured','simulated');
    ylabel('p01')
    %ylabel('lambda')
    grid on;
    title('state transition probabilities')
  subplot(212)
    handy = plot([1:L-1],p10,'r',[1:L-1],p10_test);
    %handy = plot([1:L-1],p10.*(1+p01./p10),'r',[1:L-1],p10_test.*(1+p01_test./p10_test));
    axe=axis; axis([0,L,axe(3:4)]);
    legend(handy,'measured','simulated');
    grid on;
    ylabel('p10')
    %ylabel('gamma')
    xlabel('lag [baud]');
 end;
 if plot_states,
  figure(plot_states);
  Mplot = min(50,Msim);			% number of realizations to plot
  subplot(211)
    colormap('Gray')
    imagesc(tapoff(1:Mplot,:)) 
    ylabel('realization')
    title('coefficient states: measured')
  subplot(212)
    colormap('Gray')
    imagesc(soff(1:Mplot,:)) 
    ylabel('realization')
    title('coefficient states: simulated')
    xlabel('lag [baud]');
 end;
 %p01avg = mean(p01(50:end)); p10avg = mean(p10(50:end));
 %Pavg = [1-p01avg,p01avg;p10avg,1-p10avg];
 %tapon_sim = zeros(M,L);
 %for m=1:M, tapon_sim(m,:) = 2-hmmgenerate(L,Pavg,eye(2)); end;
 %tapoff_sim = 1-tapon_sim;

 % analyze on-tap coefficient correlation
 %
 %%this is one way to do it...
 %lag_max = 2;
 %corr_h = NaN*ones(lag_max,L-1);
 %for lag = 1:lag_max,
 %  for l=1:L-lag,
 %    taponon = intersect( find(tapon(:,l)>0), find(tapon(:,l+lag)>0) ); 
 %    corr_h(lag,l) = corr(h(taponon,l),h(taponon,l+lag));
 %  end;
 %end;
 %
 %%this tests another way to do it...
 %jim = 0.2;
 %h_test(:,1) = randn(M,2)*[1;1i]/sqrt(2);
 %for l=2:L,
 %  h_test(:,l) = jim*h_test(:,l-1) + sqrt(1-jim^2)*randn(M,2)*[1;1i]/sqrt(2);
 %end;
 %corr_test = corr(h_test.*tapon); % notice restriction to "on" taps
 %handy=plot([1:L],diag(real(corr_test),1));
 %hold on;
 %  handy=[handy; plot([1:L-1],diag(imag(corr_test),1),'r')];
 %  handy=[handy; plot([1:L-2],diag(real(corr_test),2),'--')];
 %  handy=[handy; plot([1:L-2],diag(imag(corr_test),2),'r--')];
 %hold off;
 %
 corr_h = corr(h.*tapon);	% restriction to "on" taps doesn't matter?
 if plot_corr,
   figure(plot_corr); clf;
   handy=plot([1:L-1],diag(real(corr_h),1));
   hold on;
    handy=[handy; plot([1:L-1],diag(imag(corr_h),1),'r')];
    handy=[handy; plot([1:L-2],diag(real(corr_h),2),'--')];
    handy=[handy; plot([1:L-2],diag(imag(corr_h),2),'r--')];
   hold off;
   legend(handy,'lag 1 (real)','lag 1 (imag)','lag 2 (real)','lag 2 (imag)');
   axe=axis; axis([0,L,axe(3:4)]);
   grid on;
   title('on-coefficient correlation')
   ylabel('correlation coefficient')
   xlabel('lag [baud]');
 end;

 if print_plots,
   if plot_em, figure(plot_em); print('-depsc2','figs/plot_em.eps'); end;
   if plot_taphist, figure(plot_taphist); print('-depsc2','figs/plot_taphist.eps'); end;
   if plot_stats, figure(plot_stats); print('-depsc2','figs/plot_stats.eps'); end;
   if plot_corr, figure(plot_corr); print('-depsc2','figs/plot_corr.eps'); end;
   if plot_mc, figure(plot_mc); print('-depsc2','figs/plot_mc.eps'); end;
   if plot_states, figure(plot_states); print('-depsc2','figs/plot_states.eps'); end;
   if plot_chan_lin, figure(plot_chan_lin); print('-depsc2','figs/plot_chan_lin.eps'); end;
   if plot_chan_dB, figure(plot_chan_dB); print('-depsc2','figs/plot_chan_dB.eps'); end;
 end;

 return

 % used to generate the "gm2.eps" plot for the JSTSP submission
 figure(9)
  subplot(311)
    %plot([0:L-1],1-lam_gm,[Lpre,Lpre],[0,1],'r--');
    plot([0:L-1],p1,[Lpre,Lpre],[0,1],'r--');
    axis([0,L,0,1]);
    grid on;
    ylabel('lambda')
  subplot(312)
    plot([1:L-1],p01,[Lpre,Lpre],[0,1],'r--');
    axis([0,L,0,1]);
    grid on;
    ylabel('p01')
  subplot(313)
    plot([1:L-1],p10,[Lpre,Lpre],[0,1],'r--');
    axis([0,L,0,1]);
    grid on;
    ylabel('p10')
    xlabel('lag [baud]');
