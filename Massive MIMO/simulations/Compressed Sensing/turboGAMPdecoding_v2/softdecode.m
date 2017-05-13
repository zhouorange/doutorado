function [LLR_cc_decpost,success,LLR_bb_from_del,LLR_bb_eqpost] = softdecode(...
		y,z_hat,mu_z,mu_v,p_s_to_y,...
		LLR_bb_to_del,s_hyp,s_hyp2,k_one,k_zero,...
		intlv,HH,LLRmax_ldpc,iter_ldpc,LLRmax_eq,...
		T,M,Mc,Mb,Nd,n_data_,mm_data_,mm_train_) 

  % initialize
  p_s_from_y = cell(1,T);
  %s_hard = cell(1,T);
  LLR_cc_from_del = NaN*ones(Mc,1);
  LLR_bb_eqpost = NaN*ones(Mb,1);

  % equalizer outputs -> symbol probabilities
  for t=1:T,
    p_s_from_y{t} = exp( -(abs(y{t}*ones(1,2^M)-z_hat{t}*s_hyp).^2 ...
                   )./(mu_z{t}*s_hyp2+mu_v*ones(1,2^M))...
                )./(mu_z{t}*s_hyp2+mu_v*ones(1,2^M));     % ~likelihoods
    p_s_from_y{t} = p_s_from_y{t}./(sum(p_s_from_y{t},2)*ones(1,2^M));
    %[dum,k_hard] = max(p_s_from_y{t}.');               % optional!
    %s_hard{t} = s_hyp(k_hard).';                       % optional!
  end;%t

  % symbol probabilities -> bit probabilities
  for t=1:T,
    p_s_tofrom_y = p_s_to_y{t}.*p_s_from_y{t};
    p_b_one = NaN*ones(M*Nd,1);                         % proportional to Pr{b=1}
    p_b_zero = NaN*ones(M*Nd,1);                        % proportional to Pr{b=0}    
    for m=1:M,
      p_b_one([0:M:M*Nd-1]+m) = sum(p_s_tofrom_y(n_data_,k_one(m,:)),2);
      p_b_zero([0:M:M*Nd-1]+m) = sum(p_s_tofrom_y(n_data_,k_zero(m,:)),2);
    end;
    LLR_bb_eqpost((t-1)*M*Nd+[1:Nd*M]) = log((p_b_one+eps)./(p_b_zero+eps));
  end;%t
  LLR_bb_from_del = LLR_bb_eqpost - LLR_bb_to_del;
  LLR_bb_from_del = min(max(LLR_bb_from_del,-LLRmax_eq),LLRmax_eq);
  % I TRIED LIMITING THE LLR_bb_eqpost BUT THIS DIDNT SEEM TO WORK WELL!
  if ~isreal(LLR_bb_from_del), error('LLR_bb_from_del is complex'); end;

  % refresh training bits
  LLR_bb_from_del(mm_train_) = LLR_bb_to_del(mm_train_);

  % LDPC decoding
  LLR_cc_from_del(intlv) = LLR_bb_from_del(mm_data_);   % deinterleave coded bits
  p_cc_from_del = 1./(1+exp(-LLR_cc_from_del));
  [LLR_cc_decpost,success] = ldpc_decode1(0,1-p_cc_from_del,p_cc_from_del,HH,...
        LLRmax_ldpc,iter_ldpc);                         % decode
