% BP-based sparsity decoding for 2-state markov process, 
% Version 4.0, 15-Jun-2011, coded by Phil Schniter, OSU.
%
% (this differs from mc_decode.m in that it uses likelihoods rather than probabilities)
%
 
function [Le, Lp] = mc_decode2(llr_a, lam_s, p10, p01)
%nargin=0;

 %%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % parameter setup
 if nargin==0,
   N = 512; 			% signal length
   sig2 = 0.2;			% CAWGN noise variance
   plot_hist = 1;
   plot_pattern = 2;
   plot_scatter = 3;
   p01 = 0.1;			% switching probability p01
   p10 = 0.025;			% switching probability p10
     P=[1-p01,p01;p10,1-p10];
     [Evec,Eval] = eig(P');	% steady-state statistics
     evec = Evec(:,find(diag(Eval)==1));
     evec_s = evec/sum(evec);	% steady-state probabilities
   lam_s = evec_s(1)*ones(N,1);	% vectorize
   p01 = p01*ones(N-1,1);	% vectorize
   p10 = p10*ones(N-1,1);	% vectorize
 else%nargin~=0
   N = length(llr_a);
   llr_a = llr_a(:);
   p10 = p10(:);
   p01 = p01(:);
   lam_s = lam_s(:);
   plot_hist = 0;
   plot_pattern = 0;
   plot_scatter = 0;
 end;%nargin

 % markov setup
 if (length(p10)~=N-1)||(length(p01)~=N-1), 
  error('p01 and p10 must be vectors of length N-1!'); 
 end;

 % input message setup
 if nargin==0,
   % generate binary markov process 
   s_true = NaN*ones(N,1);
   s_true(1) = (rand(1)<lam_s(1));	% initialize
   for n=1:N-1,				% simulate markov chain
     s_true(n+1) = (rand(1) < [s_true(n),1-s_true(n)]*[1-p01(n);p10(n)]);
   end;
   S_true = find(s_true==1);
   S_true_not = find(s_true==0);
   if plot_pattern,
     figure(plot_pattern); clf;
     ss = NaN*ones(ceil(sqrt(N))); ss(1:N) = s_true; imagesc(ss)
     title('pattern realization')
   end;

   % generate unscaled observations in CAWGN
   y = s_true + randn(N,2)*[1;1i]*sqrt(sig2/2);
   thresh = (1+sig2*log(1./lam_s-1))/2;	% MAP detection threshold assuming iid source
   s_hat_thresh = (real(y)>thresh);	% MAP bit estimates assuming iid source
   if plot_scatter,
     figure(plot_scatter); clf;
     plot(real(y(S_true_not)),imag(y(S_true_not)),'b.',...
     	real(y(S_true)),imag(y(S_true)),'r.',...
   	[0,1],[0,0],'k+'); axe = axis;
     hold on; plot([1,1]*thresh(1),axe(3:4),'k--'); hold off;% valid only for constant thresh!
     title('observations'); xlabel('real'); ylabel('imag');
   end;
   % generate input message
   mu_f_to_s = exp(-abs(y*[1,1]-[ones(N,1),zeros(N,1)]).^2/sig2)/pi/sig2;	% unscaled
   lr_a = mu_f_to_s(:,1)./mu_f_to_s(:,2);

 else%nargin~=0
   % generate input message
   lr_a = exp(llr_a);

%------------------------------------
% only used for old way of doing things
mu_f_to_s = zeros(N,2);
mu_f_to_s(:,1) = 1./(1+exp(-llr_a));
mu_f_to_s(:,2) = 1-mu_f_to_s(:,1);						% scaled
%------------------------------------

 end;%nargin

 %%% BP INFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % forward pass
 lr_g_to_s_forward = zeros(N,1);
 lr_g_to_s_forward(1) = lam_s(1)/(1-lam_s(1)+eps);
 for n=2:N,
   lr_g_to_s_forward(n) = ( lr_g_to_s_forward(n-1)*lr_a(n-1)*(1-p01(n-1)) + p10(n-1) ...
   	)/( lr_g_to_s_forward(n-1)*lr_a(n-1)*p01(n-1) + (1-p10(n-1)) + eps );
 end;

 % backward pass 
 lr_g_to_s_backward = zeros(N-1,1);
 lr_g_to_s_backward(N-1) = ( lr_a(N)*(1-p01(N-1)) + p01(N-1) ...
 	)/( lr_a(N)*p10(N-1) + (1-p10(N-1)) );
 for n=N-2:-1:1,
   lr_g_to_s_backward(n) = ( lr_g_to_s_backward(n+1)*lr_a(n+1)*(1-p01(n)) + p01(n) ...
   	)/( lr_g_to_s_backward(n+1)*lr_a(n+1)*p10(n) + (1-p10(n)) + eps );
 end;

 % output messages
 lr_s_to_f = lr_g_to_s_forward.*[lr_g_to_s_backward;1];
 Le = log(lr_s_to_f);

 % bit inference
 Lp = log(lr_a) + Le;


%------------------------------------
if 0,
 % forward pass: old way
 mu_g_to_s_forward = zeros(N,2);
 mu_g_to_s_forward(1,:) = [lam_s(1),1-lam_s(1)];
 for n=2:N,
   P=[1-p01(n-1),p01(n-1);p10(n-1),1-p10(n-1)];
   mu_gn_to_sn_forward = ( mu_g_to_s_forward(n-1,:).*mu_f_to_s(n-1,:) )*P;
   mu_g_to_s_forward(n,:) = mu_gn_to_sn_forward/sum(mu_gn_to_sn_forward);
 end;
 
 % backward pass: old way
 mu_g_to_s_backward = zeros(N-1,2);
 mu_gn_to_sn_backward = mu_f_to_s(N,:)*[1-p01(N-1),p01(N-1);p10(N-1),1-p10(N-1)].';
 mu_g_to_s_backward(N-1,:) = mu_gn_to_sn_backward/sum(mu_gn_to_sn_backward);
 for n=N-2:-1:1,
   Pt=[1-p01(n),p01(n);p10(n),1-p10(n)].';
   mu_gn_to_sn_backward = ( mu_g_to_s_backward(n+1,:).*mu_f_to_s(n+1,:) )*Pt;
   mu_g_to_s_backward(n,:) = mu_gn_to_sn_backward/sum(mu_gn_to_sn_backward);
 end;
 
 % output messages: old way
 mu_s_to_f = mu_g_to_s_forward.*[mu_g_to_s_backward;0.5,0.5];
 Le_old = log(mu_s_to_f(:,1)./mu_s_to_f(:,2));
 
 % bit inference: old way
 p_s_given_y = mu_s_to_f.*mu_f_to_s;
 %p_s_given_y = p_s_given_y./(p_s_given_y*ones(2));	% normalize (optional)
 Lp_old = log(p_s_given_y(:,1)./p_s_given_y(:,2));
 
 % check stuff for debugging
 sum(abs(Le(2:end)-Le_old(2:end)))
 sum(abs(Lp(2:end)-Lp_old(2:end)))
 %[Le,Le_old,Lp,Lp_old]
end;
%------------------------------------

 %%% REPORT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if nargin==0,
   % error rate
   s_hat = (Lp>0);
   Pe_bp = mean(s_hat~=s_true) 
   Pe_thresh = mean(s_hat_thresh~=s_true) 
   
   if plot_hist, 
     figure(plot_hist); clf;
     %[dum,llr_bins] = hist(Le.*(2*s_true-1),50);
     %[dum,llr_bins] = hist(Le,51);
     llr_bins = max(abs(Le))*linspace(-1.1,1.1,100);
     hts_on = hist(Le(S_true),llr_bins)/length(S_true);
     hts_off = hist(Le(S_true_not),llr_bins)/length(S_true_not);
     bar(llr_bins,[hts_off;hts_on]','hist'); axe=axis;
     hold on; plot([0,0],axe(3:4),'k--'); hold off;
     title('extrinsic LLR pdf')
   end;
 end;

