%%% CompChaEst.m
%%%
%%% Simulates data transmission over a wireless doubly selective/dispersive
%%% MIMO channel using an OFDM scheme, where channel estimation is done using
%%% algorithms from compressed sensing, modified compressed sensing,
%%% distributed compressed sensing and group-sparse recovery algorithms
%%%
%%%
%%% COPYRIGHT : (c) Daniel Eiwen, Georg Tauböck, 2009-2012
%%%
%%% For details see
%%% D. Eiwen, "Compressive Channel Estimation - Compressed Sensing Methods
%%% for Estimating Doubly Selective Channels in Multicarrier Systems",
%%% PhD-thesis, University of Vienna, Vienna, May 2012
%%%
%%% The toolboxes cvx, spgl and IlmProp have to be installed for this file to work (see Readme).

clear all
format short e
rand('state', sum(100*clock));


%%% Simulation parameters

Ges_iter=1;     %% number of iterations
NOISE_LEVEL=[0.1, 0.02, 0.01, 0.005, 0.001];  %% relative level of AWGN

Track=1;        %% set to 1 if tracking should be performed, and 0 otherwise
steps=10;       %% number of time steps per iteration for tracking
suppsize=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];  %% sizes of the parts 
                %%% of the support that are kept (0 meaning that no inform-
                %%% ation, and 1 meaning that the whole support is kept)
lsupp=length(suppsize);
if Track==0 
    steps=1; 
end

%%% Recovery algorithms (set 1 if the algorithm should be used)

LS=1;               %% Least squares channel estimation
%%% Algorithms for channel-per-channel estimation
%%% Note that in the case of group size (1,1) (which corresponds to
%%% gr_i=gr_m=0) these algorithms coincide with their conventional
%%% counterparts.
GBPDN_DFT=1;        %% Group-BPDN using DFT basis
GBPDN_OPT=1;        %% Group-BPDN using optimized basis
GCOSAMP_DFT=1;      %% Group-CoSaMP using DFT basis
GCOSAMP_OPT=1;      %% Group-CoSaMP using optimized basis
GOMP_DFT=1;         %% Group-OMP using DFT basis
GOMP_OPT=1;         %% Group-OMP using optimized basis
%%% Algorithms for multichannel estimation
%%% Note that these algorithms coincide with their conventional
%%% counterparts (GBPDN, GCOSAMP and GOMP) in the SISO case.
MGBPDN_DFT=1;        %% Group-BPDN using DFT basis
MGBPDN_OPT=1;        %% Group-BPDN using optimized basis
MGCOSAMP_DFT=1;      %% Group-CoSaMP using DFT basis
MGCOSAMP_OPT=1;      %% Group-CoSaMP using optimized basis
GDCSSOMP_DFT=1;      %% Group-DCSSOMP using DFT basis
GDCSSOMP_OPT=1;      %% Group-DCSSOMP using optimized basis
%%% Algorithms for channel tracking
%%% Note that in the case that no support information is given (i.e. if
%%% suppsize=0), MODGDCSSOMP coincides with conventional G-DCSSOMP.
MODGDCSSOMP_DFT=1;   %% Modified Group-DCSSOMP using DFT basis
MODGDCSSOMP_OPT=1;   %% Modified Group-DCSSOMP using optimized basis

%%% Group sizes (exponent of 2)
gr_i=[0,1,2];          %% group sizes in i-direction (exponent of 2)
gr_m=[0,1];            %% group sizes in m-direction (exponent of 2)
ginum=length(gr_i);    %% used for counting
gmnum=length(gr_m);    %% used for counting

%%% Parameters for the reconstruction algorithms

%%% Channel-per-channel estimation
it_gbpdn_dft=10^(-2.6);                  %% G-BPDN parameter using DFT basis
it_gbpdn_opt=10^(-2.6);                  %% G-BPDN parameter using optimized basis
it_gcosamp_dft=[100,75;65,38;42,22];     %% sparsity estimate for G-COSAMP using DFT basis
it_gcosamp_opt=[80,70;50,28;34,20];      %% sparsity estimate for G-COSAMP using optimized basis
it_gomp_dft=[120,100;80,44;38,21];       %% number of G-OMP iterations using DFT basis
it_gomp_opt=[100,85;55,31;42,24];        %% number of G-OMP iterations using optimized basis
%%% Multichannel estimation
it_mgbpdn_dft=10^(-2.6);                 %% G-BPDN parameter using DFT basis
it_mgbpdn_opt=10^(-2.6);                 %% G-BPDN parameter using optimized basis
it_mgcosamp_dft=[130,90;75,40;45,23];    %% sparsity estimate for G-COSAMP using DFT basis
it_mgcosamp_opt=[80,65;45,25;28,17];     %% sparsity estimate for G-COSAMP using optimized basis
it_gdcssomp_dft=[130,90;80,52;47,27];    %% number of G-DCSSOMP iterations using DFT basis
it_gdcssomp_opt=[110,70;60,30;28,18];    %% number of G-DCSSOMP iterations using optimized basis
%%% Channel tracking
it_modgdcssomp_dft=zeros (ginum,gmnum,lsupp);   %% number of MODGDCSSOMP iterations using DFT basis
for kk=1:9
    it_modgdcssomp_dft(:,:,kk)=[130,90;80,52;47,27];
end
it_modgdcssomp_dft(:,:,10)=[135,94;84,54;50,30];
it_modgdcssomp_dft(:,:,11)=[140,98;88,56;52,32];

it_modgdcssomp_opt=zeros (ginum,gmnum,lsupp);   %% number of MODGDCSSOMP iterations using optimized basis
for kk=1:9
    it_modgdcssomp_opt(:,:,kk)=[110,70;60,30;28,18];
end
it_modgdcssomp_opt(:,:,10)=[115,74;64,32;30,20];
it_modgdcssomp_opt(:,:,11)=[115,74;64,32;30,20];

%%% Note that the choice of these parameters for the algorithms is crucial
%%% for estimation quality, the values given are just an example for an
%%% appropriate order of magnitude. For more details see [1].


%%% System parameters

K=512;              %% number of subcarriers
L_cp=128;           %% length of cyclic prefix
L=32;               %% number of OFDM-symbols
N=K+L_cp;           %% length of one OFDM-symbol
tau_max=L_cp;       %% maximal delay of the channel
Td=N*L;             %% overall length of data transmission vector
bandwidth=5e6;      %% bandwith
f0=5e9;             %% carrier frequency

NT=2                %% number of transmit antennas
NR=2                %% number of receive antennas
Theta=NR*NT;        %% overall number of cross channels

g1=ones(1,N);       %% transmit pulse for CP-OFDM
gamma1=[zeros(1,L_cp) ones(1,K)];   %% receive pulse for CP-OFDM

delta_L=1;          %% subgrid spacing delta L
delta_K=K/(tau_max);%% subgrid spacing delta K

J=L/delta_L;
D=K/delta_K;


%%% Selecting pilots

Numb_pilots=1024;                          %% number of pilots per step 
pilots=sqrt(NT)*(1+1j)*eye(NT);            %% "diagonal" pilot matrix
% pilots=sqrt(NT*2)*orth(randc(NT));       %% random, orthonormal pilots

%%% LS channel estimation only works with the diagonal pilot setup!
%%% Note also that LS channel estimation needs a lot more pilots than compressive channel estimation!
if LS==1 && norm(pilots-diag(diag(pilots)),'fro')>0 
    disp('');
    disp('LS channel estimation only works with the diagonal pilot setup! Therefore LS estimation will not be performed.');
    tst=1;
    while tst
        rep=input('Do you want to continue? (y/n): ','s');
        if strcmp(rep,'y')
            LS=0;
            tst=0;
            else if strcmp(rep,'n')
                return;
            end
        end
    end
end

inv_pilots=inv(pilots);
pil_pos=zeros(steps*Numb_pilots,2,NT);
pilot_positions=zeros(steps*Numb_pilots,1);

for st=1:steps
    permute_vec=randperm(J*D).';	     %% here we use different pilot positions for the different time steps
    for ss=1:NT
        pp=sort(permute_vec(((ss-1)*Numb_pilots+1):ss*Numb_pilots));
        for rr=1:NR
            pilot_positions((st-1)*Numb_pilots+1:st*Numb_pilots,(ss-1)*NR+rr)=pp;
        end
        pil_pos((st-1)*Numb_pilots+1:st*Numb_pilots,2,ss)=1+floor((pp-1)./D)*delta_L;
        pil_pos((st-1)*Numb_pilots+1:st*Numb_pilots,1,ss)=1+mod((pp-1),D)*delta_K;
    end
end

Help_pilots=zeros(K,steps*L,NT);        %% just for simplification later on
for kk=1:Numb_pilots
    for st=1:steps
        for ss=1:NT
            for rr=1:NR
                Help_pilots(pil_pos((st-1)*Numb_pilots+kk,1,rr),(st-1)*L+pil_pos((st-1)*Numb_pilots+kk,2,rr),ss)=10*rr+ss;
            end
        end
    end
end
Help_pilots=sparse(Help_pilots(:));



%%% Basis optimization
disp('basis optimization');

%%% Basis optimization for channel-per-channel estimation/tracking
if GBPDN_OPT+GCOSAMP_OPT+GOMP_OPT>0
    U_SISO=zeros(J,J,ginum);
    for gi=1:ginum
        U_SISO(:,:,gi) = optimize_basis_jointgroup(gamma1,g1,J,L,K,N,0.03,-0.03,1,0,2^gr_i(gi))';
    end
%     load('U_SISO.mat');        %% alternatively, precalculate U_SISO and load it here
%                                %% in U_SISO.mat, the array U_SISO has to be
%                                %% of size J x J x ginum
    Uhat=fft(U_SISO);           %% needed for reconstruction later
    U_SISO_hat=zeros(L,J,ginum);
    U_SISO_hat(1:J/2,:,:)=Uhat(1:J/2,:,:);
    U_SISO_hat(L-J/2+1:L,:,:)=Uhat(J/2+1:J,:,:);

    normalize_vec_SISO=zeros(J*D,Theta,ginum);   %% vector of normalization values for optimized basis, i.e. the diagonal of the scaling matrices D^(q)
    for tt=1:Theta
        for gi=1:ginum
            for kk=1:D*J
                ekk=zeros(D*J,1);
                ekk(kk)=1;
                for st=1:steps
                    normalize_vec_SISO(kk,tt,gi,st)=1/(norm(A_optimized(ekk,pilot_positions((st-1)*Numb_pilots+1:st*Numb_pilots,:),U_SISO(:,:,gi),D,J,ones(D*J,Theta),tt,1)));
                end
            end
        end
    end
end

%%% Basis optimization for multichannel estimation/tracking
if MGBPDN_OPT+MGCOSAMP_OPT+GDCSSOMP_OPT+MODGDCSSOMP_OPT>0
    U_MIMO=zeros(J,J,ginum);
    for gi=1:ginum
        U_MIMO(:,:,gi) = optimize_basis_jointgroup(gamma1,g1,J,L,K,N,0.03,-0.03,Theta,3,2^gr_i(gi))';
    end
%     load('U_MIMO.mat');       %% alternatively, precalculate U_MIMO and load it here
%                               %% in U_MIMO.mat, the array U_MIMO has to be
%                               %% of size J x J x ginum
    Uhat=fft(U_MIMO);           %% needed for reconstruction later
    U_MIMO_hat=zeros(L,J,ginum);
    U_MIMO_hat(1:J/2,:,:)=Uhat(1:J/2,:,:);
    U_MIMO_hat(L-J/2+1:L,:,:)=Uhat(J/2+1:J,:,:);
    
    normalize_vec_MIMO=zeros(J*D,Theta,ginum);   %% vector of normalization values for optimized basis, i.e. the diagonal of the scaling matrices D^(q)
    for tt=1:Theta
        for gi=1:ginum
            for kk=1:D*J
                ekk=zeros(D*J,1);
                ekk(kk)=1;
                for st=1:steps
                    normalize_vec_MIMO(kk,tt,gi,st)=1/(norm(A_optimized(ekk,pilot_positions((st-1)*Numb_pilots+1:st*Numb_pilots,:),U_MIMO(:,:,gi),D,J,ones(D*J,Theta),tt,1)));
                end
            end
        end
    end
end

disp('basis optimization done');
disp('');


nl=length(NOISE_LEVEL);

%%% Initialization

MSE_all_N=0;
MSE_N=zeros(steps,1);

if LS==1, MSE_LS=zeros(nl,1); end

%%% Initialization for channel-per-channel estimation
if GBPDN_DFT==1, MSE_GBPDN_DFT=zeros(ginum,gmnum,nl); end
if GBPDN_OPT==1, MSE_GBPDN_OPT=zeros(ginum,gmnum,nl); end
if GCOSAMP_DFT==1, MSE_GCOSAMP_DFT=zeros(ginum,gmnum,nl); end
if GCOSAMP_OPT==1, MSE_GCOSAMP_OPT=zeros(ginum,gmnum,nl); end
if GOMP_DFT==1, MSE_GOMP_DFT=zeros(ginum,gmnum,nl); end
if GOMP_OPT==1, MSE_GOMP_OPT=zeros(ginum,gmnum,nl); end

%%% Initialization for multichannel estimation
if MGBPDN_DFT==1, MSE_MGBPDN_DFT=zeros(ginum,gmnum,nl); end
if MGBPDN_OPT==1, MSE_MGBPDN_OPT=zeros(ginum,gmnum,nl); end
if MGCOSAMP_DFT==1, MSE_MGCOSAMP_DFT=zeros(ginum,gmnum,nl); end
if MGCOSAMP_OPT==1, MSE_MGCOSAMP_OPT=zeros(ginum,gmnum,nl); end
if GDCSSOMP_DFT==1, MSE_GDCSSOMP_DFT=zeros(ginum,gmnum,nl); end
if GDCSSOMP_OPT==1, MSE_GDCSSOMP_OPT=zeros(ginum,gmnum,nl); end

%%% Initialization for channel tracking
if MODGDCSSOMP_DFT==1, MSE_MODGDCSSOMP_DFT=zeros(lsupp,ginum,gmnum,nl); end
if MODGDCSSOMP_OPT==1, MSE_MODGDCSSOMP_OPT=zeros(lsupp,ginum,gmnum,nl); end


%%% Main loop (over the channels)
disp('main loop');
for iter=1:Ges_iter
    
    iter
    
    alpha=0.25;             %% roll-off factor for RRC-filter
    Tds=steps*N*L;          %% overall length of data transmission vector
    LL=tau_max*(1+alpha)+1; %% delay-width for RRC-filter
    tau_start=0;            %% delay correction (-> shift to the left so that no cyclic artifacts occur)
    
    support_dft=cell(lsupp,ginum,gmnum,nl); %% support that is kept for channel tracking using the DFT basis
    support_opt=cell(lsupp,ginum,gmnum,nl); %% support that is kept for channel tracking using an optimized basis
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate channel
    disp('calculating channel');
    
    %%% Free the memory
    clear HH h spread_funct F_H F_H_sampled s r rn H_calc H_BPDN_DFT H_OMP_DFT H_SOMP_DFT H_DCSSOMP_DFT H_SOMP_OPT H_DCSSOMP_OPT
    
    H=gen_channel_cce(N,L,steps,NT,NR,bandwidth,f0,LL,iter);   %% generates a channel using the toolbox IlmProp
%     tmpl_l='channel_%2.3d';           %% alternatively, Channels can be precalculated
%     cd Channels                       %% and loaded here
%     fname_l=sprintf(tmpl_l,iter);     %% here: "channel_001.mat" ...
%     load(fname_l, 'H');
%     cd ..
    
    %%% HH ... reordering of the cross-channels from the output of IlmProp
    %%% to be compatible prepare channel H as NR x NT x Tds x LL array
    HH=zeros(Tds,LL,Theta);
    for ss=1:NT
        for rr=1:NR
            HH(:,:,(ss-1)*NR+rr)=squeeze(H(rr,ss,:,:));
        end
    end
    
    clear H;
    
    HH=permute(HH,[2,1,3]).*repmat( rcosine(LL , alpha ), [1,Tds,Theta]  ) .* repmat( exp(1j*2*pi* (-(LL-1)/2:(LL-1)/2).' *tau_start*(1+alpha)/(LL-1)),[1,Tds,Theta] );
    
    h=permute( tau_max*ifft(HH((LL+1)/2:end,:,:),tau_max) + fft(HH((LL+1)/2:-1:1,:,:),tau_max) - repmat(HH((LL+1)/2,:,:),[tau_max,1,1]) , [2,1,3]) ;
    clear HH;
    
    disp('channel calculated');
    
    spread_funct=zeros(Td,tau_max,Theta,steps); %% spreading function of the channel
    F_H=zeros(Td,tau_max,Theta,steps);          %% "Fourier coefficients" of the channel
    F_H_sampled=zeros(L,tau_max,Theta,steps);   %% sampled version
    H_calc=zeros(K,L,Theta,steps);              %% actual channel matrix
    
    for st=1:steps
        for tt=1:Theta
            spread_funct(:,:,tt,st)=fft(h((st-1)*Td+1:st*Td,:,tt)/Td);
            F_H(:,:,tt,st)=ifft(spread_funct(:,:,tt,st).*conj(ambiguity_function(gamma1,g1,tau_max,Td)));
            F_H_sampled(:,:,tt,st)=F_H(1:N:Td,:,tt,st);
            H_calc(:,:,tt,st)=(N*L/K)*(fft(F_H_sampled(:,:,tt,st).',K));   
        end
        MSE_N(st)=(sum(sum(sum((abs(H_calc(:,:,:,st))).^2))));
        if st==steps, disp('channel coefficients calculated'); end
    end
    
    
    
    %%% Data Transmission
    
    disp('start transmission');
    
    a=zeros(K,steps*L,NT);      %% transmitted symbols (QPSK)
    x=zeros(K,steps*L,NR);      %% demodulated symbols
    rn=zeros(Tds,NR);           %% noisy received signal
    r=zeros(Tds,NR);            %% noiseless received signals
    s=zeros(Tds,NT);            %% transmitted signal
    
    %%% Transmitter
    
    symbols=[1+1j,1-1j,-1+1j,-1-1j];      %% QPSK symbols
    a(:)=symbols(ceil(4*rand(K*L*NT*steps,1)));
    for rr=1:NR
        for ss=1:NT
            a(Help_pilots==10*ss+rr)=pilots(rr,ss);
        end
    end
    for ss=1:NT
        s(:,ss) = MC_modulator(a(:,:,ss),g1,N);
    end
    
    
    %%% Transmission
    
    norm_r=zeros(1,Theta);
    for tt=1:Theta
        r(:,tt) = channel_calc(tau_max,s(:,ceil(tt/NR)),h(:,:,tt));
        norm_r(tt)=norm(r(:,tt))^2;
    end
    
    for rr=1:NR
        w_hlp(rr)=sum(norm_r(rr:NR:end));
    end
    whlp=mean(w_hlp);
    
    
    %%% Inner loop over noise level (SNR)
    for reln=1:length(NOISE_LEVEL)
        rel_noise = NOISE_LEVEL(reln)
        
        %%% Receiver
        
        for rr=1:NR
            r_hlp=sum(r(:,rr:NR:end),2);
            w = randn(Tds, 1) + sqrt(-1)*randn(Tds, 1);
            w = ((sqrt(rel_noise*whlp ))  / norm(w)) * w;
            rn(:,rr) = r_hlp + w;
            x(:,:,rr)=MC_demodulator(rn(:,rr),gamma1,K,N);
        end

        %%% Channel estimation
        disp('start channel estimation');
        
        v_st=zeros(Numb_pilots,Theta,steps);     %% "vector of measurements"
        for ss=1:NT
            for rr=1:NR
                for st=1:steps
                    for kk=1:Numb_pilots
                        v_st(kk,(ss-1)*NR+rr,st)=x(pil_pos((st-1)*Numb_pilots+kk,1,ss),(st-1)*L+pil_pos((st-1)*Numb_pilots+kk,2,ss),rr);
                    end
                end
            end
        end
        
        u_dft=zeros(D*J,Theta);
        u_opt=zeros(D*J,Theta);
        
        for ll=1:nl
            for gi=1:ginum
                for gm=1:gmnum
                    for susi=1:lsupp
                        support_dft{susi,gi,gm,ll}=[];
                        support_opt{susi,gi,gm,ll}=[];
                    end
                end
            end
        end
        
        
        %%% Loop over the time steps
        for st=1:steps
            
            v=v_st(:,:,st);
            
            %%% LS channel estimation
            if LS==1
                hlpJ=zeros(J+1,NT);
                for ss=1:NT
                    for jj=1:J
                        hlpJ(jj+1,ss)=hlpJ(jj,ss)+length(find(pil_pos((st-1)*Numb_pilots+1:st*Numb_pilots,2,ss)==(jj-1)*delta_L+1));
                    end
                end
                disp('LS'); H_LS=zeros(K,L,Theta);
                for rr=1:NR
                    for ss=1:NT
                        for jj=1:J
                            H_LS(:,(jj-1)*delta_L+1,(ss-1)*NR+rr)=spline(pil_pos((st-1)*Numb_pilots+(hlpJ(jj,ss)+1:hlpJ(jj+1,ss)),1,ss),v(hlpJ(jj,ss)+1:hlpJ(jj+1,ss),(ss-1)*NR+rr)/pilots(ss,ss).',1:K).';
                        end
                        H_LS(:,:,(ss-1)*NR+rr)=spline((1:delta_L:L),H_LS(:,1:delta_L:end,(ss-1)*NR+rr),1:L);
                    end
                end
                MSE_LS(reln)=MSE_LS(reln)+(sum(sum(sum((abs(H_LS-H_calc(:,:,:,st))).^2))));
            end
        
            %%% Fourier measurement matrix
            A_dft= @(z,mat,mode) A_Fourier(z,pilot_positions((st-1)*Numb_pilots+1:st*Numb_pilots,:),D,J,mat,mode);          %% Fourier measurement matrix
            A_dft_ext= @(z,mode) A_Fourier_ext(z,pilot_positions((st-1)*Numb_pilots+1:st*Numb_pilots,:),D,J,Theta,mode);    %% Fourier measurement matrix for the extended system
            A_dft_sel=@(z,T,mat,mode)A_selectcols(z,T,A_dft,D*J,mat,mode);            %% function handle for selecting columns (necessary for gcosamp.m)
            A_dft_ext_sel=@(z,T,mode)A_selectcols_ext(z,T,A_dft_ext,Theta*D*J,mode);  %% function handle for selecting columns (necessary for gcosamp.m)
            
            for gi=1:ginum
                %%% Optimized measurement matrix for channel-per-channel estimation
                if GBPDN_OPT+GCOSAMP_OPT+GOMP_OPT>0
                    norm_vec_SISO=squeeze(normalize_vec_SISO(:,:,gi));
                    UU_SISO=squeeze(U_SISO(:,:,gi));
                    UU_SISO_hat=squeeze(U_SISO_hat(:,:,gi));
                    A_SISO= @(z,mat,mode) A_optimized(z,pilot_positions((st-1)*Numb_pilots+1:st*Numb_pilots,:),UU_SISO,D,J,norm_vec_SISO,mat,mode);
                    A_SISO_sel=@(z,T,mat,mode)A_selectcols(z,T,A_SISO,D*J,mat,mode);  %% function handle for selecting columns (necessary for gcosamp.m)
                end
                %%% Optimized measurement matrix for multichannel estimation
                if MGBPDN_OPT+MGCOSAMP_OPT+GDCSSOMP_OPT+MODGDCSSOMP_OPT>0
                    norm_vec_MIMO=squeeze(normalize_vec_MIMO(:,:,gi));  %% just for speeding up multiplications
                    UU_MIMO=squeeze(U_MIMO(:,:,gi));                    %% just for speeding up multiplications
                    UU_MIMO_hat=squeeze(U_MIMO_hat(:,:,gi));            %% just for speeding up multiplications
                    A_MIMO= @(z,mat,mode) A_optimized(z,pilot_positions((st-1)*Numb_pilots+1:st*Numb_pilots,:),UU_MIMO,D,J,norm_vec_MIMO,mat,mode);
                    A_MIMO_ext= @(z,mode) A_optimized_ext(z,pilot_positions((st-1)*Numb_pilots+1:st*Numb_pilots,:),UU_MIMO,D,J,norm_vec_MIMO,Theta,mode);
                    A_MIMO_sel=@(z,T,mode)A_selectcols(z,T,A_MIMO,D*J,mode);                    %% function handle for selecting columns (necessary for gcosamp.m)
                    A_MIMO_ext_sel=@(z,T,mode)A_selectcols_ext(z,T,A_MIMO_ext,Theta*D*J,mode);  %% function handle for selecting columns (necessary for gcosamp.m)
                end
                
                for gm=1:gmnum
                    B=vec(kron(vec2mat(1:J*D/2^(gr_i(gi)+gr_m(gm)),J/2^gr_i(gi)).',ones(2^gr_i(gi),2^gr_m(gm))).'); %% group indices for GSCS-methods
                    B_ext=repmat(B,Theta,1);	%% group indices for extended system

                    %%% Channel-per-channel estimation
                    
                    if GBPDN_DFT==1
                        disp('GBPDN_DFT'); H_GBPDN_DFT=zeros(K,L,Theta);
                        for tt=1:Theta, u_dft(:,tt) = spg_group(@(z,mode)A_dft(z,tt,mode),v(:,tt),B,it_gbpdn_dft); end
                        for kk=1:D*J, for rr=1:NR, u_dft(kk,rr:NR:end)=u_dft(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_GBPDN_DFT(:,:,tt)=sqrt(J*D/Numb_pilots)*(repmat((exp((-1j)*pi*(J/L)*((0:(L-1))'))),1,K).*(ifft(((L/sqrt(J*D))*fft((vec2mat((u_dft(:,tt)),D)).',K)).',L))).'; end
                        MSE_GBPDN_DFT(gi,gm,reln)=MSE_GBPDN_DFT(gi,gm,reln)+(sum(sum(sum((abs(H_GBPDN_DFT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if GCOSAMP_DFT==1
                        disp('GCOSAMP_DFT'), H_GCOSAMP_DFT=zeros(K,L,Theta);
                        for tt=1:Theta, u_dft(:,tt) = gcosamp(v(:,tt),@(z,mode)A_dft(z,tt,mode),@(z,T,mode)A_dft_sel(z,T,tt,mode),B,it_gcosamp_dft(gi,gm),16); end
                        for kk=1:D*J, for rr=1:NR, u_dft(kk,rr:NR:end)=u_dft(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_GCOSAMP_DFT(:,:,tt)=sqrt(J*D/Numb_pilots)*(repmat((exp((-1j)*pi*(J/L)*((0:(L-1))'))),1,K).*(ifft(((L/sqrt(J*D))*fft((vec2mat((u_dft(:,tt)),D)).',K)).',L))).'; end
                        MSE_GCOSAMP_DFT(gi,gm,reln)=MSE_GCOSAMP_DFT(gi,gm,reln)+(sum(sum(sum((abs(H_GCOSAMP_DFT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if GOMP_DFT==1
                        disp('GOMP_DFT'); H_GOMP_DFT=zeros(K,L,Theta);
                        for tt=1:Theta, u_dft(:,tt) = modgdcssomp(v(:,tt),@(z,mat,mode)A_dft(z,tt,mode),B,it_gomp_dft(gi,gm),[],0.01); end
                        %%% note that GOMP is a special case of MODGDCSSOMP for just one channel and no support information
                        for kk=1:D*J, for rr=1:NR, u_dft(kk,rr:NR:end)=u_dft(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_GOMP_DFT(:,:,tt)=sqrt(J*D/Numb_pilots)*(repmat((exp((-1j)*pi*(J/L)*((0:(L-1))'))),1,K).*(ifft(((L/sqrt(J*D))*fft((vec2mat((u_dft(:,tt)),D)).',K)).',L))).'; end
                        MSE_GOMP_DFT(gi,gm,reln)=MSE_GOMP_DFT(gi,gm,reln)+(sum(sum(sum((abs(H_GOMP_DFT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if GBPDN_OPT==1
                        disp('GBPDN_OPT'); H_GBPDN_OPT=zeros(K,L,Theta);
                        for tt=1:Theta, u_opt(:,tt) = spg_group(@(z,mode)A_SISO(z,tt,mode),v(:,tt),B,it_gbpdn_opt).*norm_vec_SISO(:,tt); end
                        for kk=1:D*J, for rr=1:NR, u_opt(kk,rr:NR:end)=u_opt(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_GBPDN_OPT(:,:,tt)=(L/J*ifft(UU_SISO_hat)*(((1/sqrt(D))*fft((vec2mat((u_opt(:,tt)),D)).',K)).')).'; end
                        MSE_GBPDN_OPT(gi,gm,reln)=MSE_GBPDN_OPT(gi,gm,reln)+(sum(sum(sum((abs(H_GBPDN_OPT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if GCOSAMP_OPT==1
                        disp('GCOSAMP_OPT'); H_GCOSAMP_OPT=zeros(K,L,Theta);
                        for tt=1:Theta, u_opt(:,tt) = gcosamp(v(:,tt),@(z,mode)A_SISO(z,tt,mode),@(z,T,mode)A_SISO_sel(z,T,tt,mode),B,it_gcosamp_opt(gi,gm),16).*norm_vec_SISO(:,tt); end
                        for kk=1:D*J, for rr=1:NR, u_opt(kk,rr:NR:end)=u_opt(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_GCOSAMP_OPT(:,:,tt)=(L/J*ifft(UU_SISO_hat)*(((1/sqrt(D))*fft((vec2mat((u_opt(:,tt)),D)).',K)).')).'; end
                        MSE_GCOSAMP_OPT(gi,gm,reln)=MSE_GCOSAMP_OPT(gi,gm,reln)+(sum(sum(sum((abs(H_GCOSAMP_OPT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if GOMP_OPT==1
                        disp('GOMP_OPT'), H_GOMP_OPT=zeros(K,L,Theta);
                        for tt=1:Theta, u_opt(:,tt) = modgdcssomp(v(:,tt),@(z,mat,mode)A_SISO(z,tt,mode),B,it_gomp_opt(gi,gm),[],0.01).*norm_vec_SISO(:,tt); end
                        for kk=1:D*J, for rr=1:NR, u_opt(kk,rr:NR:end)=u_opt(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_GOMP_OPT(:,:,tt)=(L/J*ifft(UU_SISO_hat)*(((1/sqrt(D))*fft((vec2mat((u_opt(:,tt)),D)).',K)).')).'; end
                        MSE_GOMP_OPT(gi,gm,reln)=MSE_GOMP_OPT(gi,gm,reln)+(sum(sum(sum((abs(H_GOMP_OPT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    %%% Multichannel estimation
                    
                    if MGBPDN_DFT==1
                        disp('MGBPDN_DFT'), H_MGBPDN_DFT=zeros(K,L,Theta);
                        u_dft(:) = spg_group(A_dft_ext,v(:),B_ext,it_mgbpdn_dft);
                        for kk=1:D*J, for rr=1:NR, u_dft(kk,rr:NR:end)=u_dft(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_MGBPDN_DFT(:,:,tt)=sqrt(J*D/Numb_pilots)*(repmat((exp((-1j)*pi*(J/L)*((0:(L-1))'))),1,K).*(ifft(((L/sqrt(J*D))*fft((vec2mat((u_dft(:,tt)),D)).',K)).',L))).'; end
                        MSE_MGBPDN_DFT(gi,gm,reln)=MSE_MGBPDN_DFT(gi,gm,reln)+(sum(sum(sum((abs(H_MGBPDN_DFT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if MGCOSAMP_DFT==1
                        disp('MGCOSAMP_DFT'); H_MGCOSAMP_DFT=zeros(K,L,Theta);
                        u_dft(:) = gcosamp(v(:),A_dft_ext,A_dft_ext_sel,B_ext,it_mgcosamp_dft(gi,gm),16);
                        for kk=1:D*J, for rr=1:NR, u_dft(kk,rr:NR:end)=u_dft(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_MGCOSAMP_DFT(:,:,tt)=sqrt(J*D/Numb_pilots)*(repmat((exp((-1j)*pi*(J/L)*((0:(L-1))'))),1,K).*(ifft(((L/sqrt(J*D))*fft((vec2mat((u_dft(:,tt)),D)).',K)).',L))).'; end
                        MSE_MGCOSAMP_DFT(gi,gm,reln)=MSE_MGCOSAMP_DFT(gi,gm,reln)+(sum(sum(sum((abs(H_MGCOSAMP_DFT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if GDCSSOMP_DFT==1
                        disp('GDCSSOMP_DFT'), H_GDCSSOMP_DFT=zeros(K,L,Theta);
                        u_dft=modgdcssomp(v,A_dft,B,it_gdcssomp_dft(gi,gm),[],0.01);
                        for kk=1:D*J, for rr=1:NR, u_dft(kk,rr:NR:end)=u_dft(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_GDCSSOMP_DFT(:,:,tt)=sqrt(J*D/Numb_pilots)*(repmat((exp((-1j)*pi*(J/L)*((0:(L-1))'))),1,K).*(ifft(((L/sqrt(J*D))*fft((vec2mat((u_dft(:,tt)),D)).',K)).',L))).'; end
                        MSE_GDCSSOMP_DFT(gi,gm,reln)=MSE_GDCSSOMP_DFT(gi,gm,reln)+(sum(sum(sum((abs(H_GDCSSOMP_DFT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if MGBPDN_OPT==1
                        disp('MGBPDN_OPT'); H_MGBPDN_OPT=zeros(K,L,Theta);
                        u_opt(:) = spg_group(A_MIMO_ext,v(:),B_ext,it_mgbpdn_opt);
                        u_opt=u_opt.*norm_vec_MIMO;
                        for kk=1:D*J, for rr=1:NR, u_opt(kk,rr:NR:end)=u_opt(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_MGBPDN_OPT(:,:,tt)=(L/J*ifft(UU_MIMO_hat)*(((1/sqrt(D))*fft((vec2mat((u_opt(:,tt)),D)).',K)).')).'; end
                        MSE_MGBPDN_OPT(gi,gm,reln)=MSE_MGBPDN_OPT(gi,gm,reln)+(sum(sum(sum((abs(H_MGBPDN_OPT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if MGCOSAMP_OPT==1
                        disp('MGCOSAMP_OPT'); H_MGCOSAMP_OPT=zeros(K,L,Theta);
                        u_opt(:) = gcosamp(v(:),A_MIMO_ext,A_MIMO_ext_sel,B_ext,it_mgcosamp_opt(gi,gm),16);
                        u_opt=u_opt.*norm_vec_MIMO;
                        for kk=1:D*J, for rr=1:NR, u_opt(kk,rr:NR:end)=u_opt(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_MGCOSAMP_OPT(:,:,tt)=(L/J*ifft(UU_MIMO_hat)*(((1/sqrt(D))*fft((vec2mat((u_opt(:,tt)),D)).',K)).')).'; end
                        MSE_MGCOSAMP_OPT(gi,gm,reln)=MSE_MGCOSAMP_OPT(gi,gm,reln)+(sum(sum(sum((abs(H_MGCOSAMP_OPT-H_calc(:,:,:,st))).^2))));
                    end
                    
                    if GDCSSOMP_OPT==1
                        disp('GDCSSOMP_OPT'); H_GDCSSOMP_OPT=zeros(K,L,Theta);
                        u_opt=modgdcssomp(v,A_MIMO,B,it_gdcssomp_opt(gi,gm),[],0.01).*norm_vec_MIMO;
                        for kk=1:D*J, for rr=1:NR, u_opt(kk,rr:NR:end)=u_opt(kk,rr:NR:end)*inv_pilots.'; end, end
                        for tt=1:Theta, H_GDCSSOMP_OPT(:,:,tt)=(L/J*ifft(UU_MIMO_hat)*(((1/sqrt(D))*fft((vec2mat((u_opt(:,tt)),D)).',K)).')).'; end
                        MSE_GDCSSOMP_OPT(gi,gm,reln)=MSE_GDCSSOMP_OPT(gi,gm,reln)+(sum(sum(sum((abs(H_GDCSSOMP_OPT-H_calc(:,:,:,st))).^2))));
                    end

                    %%% Channel tracking

                    for susi=1:lsupp
                        
                        if MODGDCSSOMP_DFT==1;
                            disp('MODGDCSSOMP_DFT')
                            H_MODGDCSSOMP_DFT=zeros(K,L,Theta);
                            u_dft = modgdcssomp(v,A_dft,B,it_modgdcssomp_dft(gi,gm,susi),support_dft{susi,gi,gm,reln},0.01);
                            uhlp=zeros(J*D/2^(gr_i(gi)+gr_m(gm)),1);
                            for kk=1:J*D/2^(gr_i(gi)+gr_m(gm)), uhlp(kk)=norm(u_dft(B==kk,:),'fro'); end
                            [us,uidx]=sort(abs(uhlp),'descend');
                            ls=length(find(uhlp));
                            for kk=1:D*J, for rr=1:NR, u_dft(kk,rr:NR:end)=u_dft(kk,rr:NR:end)*inv_pilots.'; end, end
                            for tt=1:Theta, H_MODGDCSSOMP_DFT(:,:,tt)=sqrt(J*D/Numb_pilots)*(repmat((exp((-1j)*pi*(J/L)*((0:(L-1))'))),1,K).*(ifft(((L/sqrt(J*D))*fft((vec2mat((u_dft(:,tt)),D)).',K)).',L))).'; end
                            MSE_MODGDCSSOMP_DFT(susi,gi,gm,reln)=MSE_MODGDCSSOMP_DFT(susi,gi,gm,reln)+(sum(sum(sum((abs(H_MODGDCSSOMP_DFT-H_calc(:,:,:,st))).^2))));
                            supp=[]; for kk=1:floor(suppsize(susi)*ls), supp=[supp; uidx(kk)]; end
                            support_dft{susi,gi,gm,reln}=supp;
                        end
                        
                        
                        if MODGDCSSOMP_OPT==1
                            disp('MODGDCSOMP_OPT');
                            H_MODGDCSSOMP_OPT=zeros(K,L,Theta);
                            u_opt = modgdcssomp(v,A_MIMO,B,it_modgdcssomp_opt(gi,gm,susi),support_opt{susi,gi,gm,reln},0.01);
                            uhlp=zeros(J*D/2^(gr_i(gi)+gr_m(gm)),1);
                            for kk=1:J*D/2^(gr_i(gi)+gr_m(gm)), uhlp(kk)=norm(u_opt(B==kk,:),'fro'); end
                            [us,uidx]=sort(abs(uhlp),'descend');
                            ls=length(find(uhlp));
                            u_opt=u_opt.*norm_vec_MIMO;
                            for kk=1:D*J, for rr=1:NR, u_opt(kk,rr:NR:end)=u_opt(kk,rr:NR:end)*inv_pilots.'; end, end
                            for tt=1:Theta, H_MODGDCSSOMP_OPT(:,:,tt)=(L/J*ifft(UU_MIMO_hat)*(((1/sqrt(D))*fft((vec2mat((u_opt(:,tt)),D)).',K)).')).'; end
                            MSE_MODGDCSSOMP_OPT(susi,gi,gm,reln)=MSE_MODGDCSSOMP_OPT(susi,gi,gm,reln)+(sum(sum(sum((abs(H_MODGDCSSOMP_OPT-H_calc(:,:,:,st))).^2))));
                            supp=[]; for kk=1:floor(suppsize(susi)*ls), supp=[supp; uidx(kk)]; end
                            support_opt{susi,gi,gm,reln}=supp;
                        end
                        
                    end  %% susi
                end  %% gm
            end  %% gi
            
            if reln==1, MSE_all_N=MSE_all_N+MSE_N(st); end
            
            %%% Save results
            
            sp.Ges_iter=Ges_iter;
            sp.iter=iter;
            sp.L=L;
            sp.K=K;
            sp.L_cp=L_cp;
            sp.N=N;
            sp.tau_max=tau_max;
            sp.g1=g1;
            sp.gamma1=gamma1;
            sp.NOISE_LEVEL=NOISE_LEVEL;
            sp.delta_L=delta_L;
            sp.delta_K=delta_K;
            sp.Numb_pilots=Numb_pilots;
            sp.NT=NT;
            sp.NR=NR;
            sp.gr_i=gr_i;
            sp.gr_m=gr_m;
            sp.suppsize=suppsize;
            sp.Track=Track;
            sp.steps=steps;

            %%% Information about recovery algorithms
            sp.LS=LS;
            sp.GBPDN_DFT=GBPDN_DFT;
            sp.GBPDN_OPT=GBPDN_OPT;
            sp.GCOSAMP_DFT=GCOSAMP_DFT;
            sp.GCOSAMP_OPT=GCOSAMP_OPT;
            sp.GOMP_DFT=GOMP_DFT;
            sp.GOMP_OPT=GOMP_OPT;
            sp.MGBPDN_DFT=MGBPDN_DFT;
            sp.MGBPDN_OPT=MGBPDN_OPT;
            sp.MGCOSAMP_DFT=MGCOSAMP_DFT;
            sp.MGCOSAMP_OPT=MGCOSAMP_OPT;
            sp.GDCSSOMP_DFT=GDCSSOMP_DFT;
            sp.GDCSSOMP_OPT=GDCSSOMP_OPT;
            sp.MODGDCSSOMP_DFT=MODGDCSSOMP_DFT;
            sp.MODGDCSSOMP_OPT=MODGDCSSOMP_OPT;
            
            sp.MSE_all_N=MSE_all_N;
            %%% Results for channel-per-channel estimation
            if LS==1, sp.MSE_LS=10*log10(MSE_LS/MSE_all_N); end
            if GBPDN_DFT==1, sp.MSE_GBPDN_DFT=10*log10(MSE_GBPDN_DFT/MSE_all_N); sp.it_gbpdn_dft=it_gbpdn_dft; end
            if GBPDN_OPT==1, sp.MSE_GBPDN_OPT=10*log10(MSE_GBPDN_OPT/MSE_all_N); sp.it_gbpdn_opt=it_gbpdn_opt; end
            if GCOSAMP_DFT==1, sp.MSE_GCOSAMP_DFT=10*log10(MSE_GCOSAMP_DFT/MSE_all_N); sp.it_gcosamp_dft=it_gcosamp_dft; end
            if GCOSAMP_OPT==1, sp.MSE_GCOSAMP_OPT=10*log10(MSE_GCOSAMP_OPT/MSE_all_N); sp.it_gcosamp_opt=it_gcosamp_opt; end
            if GOMP_DFT==1, sp.MSE_GOMP_DFT=10*log10(MSE_GOMP_DFT/MSE_all_N); sp.it_gomp_dft=it_gomp_dft; end
            if GOMP_OPT==1, sp.MSE_GOMP_OPT=10*log10(MSE_GOMP_OPT/MSE_all_N); sp.it_gomp_opt=it_gomp_opt; end
            %%% Results for multichannel estimation
            if MGBPDN_DFT==1, sp.MSE_MGBPDN_DFT=10*log10(MSE_MGBPDN_DFT/MSE_all_N); sp.it_mgbpdn_dft=it_mgbpdn_dft; end
            if MGBPDN_OPT==1, sp.MSE_MGBPDN_OPT=10*log10(MSE_MGBPDN_OPT/MSE_all_N); sp.it_mgbpdn_opt=it_mgbpdn_opt; end
            if MGCOSAMP_DFT==1, sp.MSE_MGCOSAMP_DFT=10*log10(MSE_MGCOSAMP_DFT/MSE_all_N); sp.it_mgcosamp_dft=it_mgcosamp_dft; end 
            if MGCOSAMP_OPT==1, sp.MSE_MGCOSAMP_OPT=10*log10(MSE_MGCOSAMP_OPT/MSE_all_N); sp.it_mgcosamp_opt=it_mgcosamp_opt; end
            if GDCSSOMP_DFT==1, sp.MSE_GDCSSOMP_DFT=10*log10(MSE_GDCSSOMP_DFT/MSE_all_N); sp.it_gdcssomp_dft=it_gdcssomp_dft; end
            if GDCSSOMP_OPT==1, sp.MSE_GDCSSOMP_OPT=10*log10(MSE_GDCSSOMP_OPT/MSE_all_N); sp.it_gdcssomp_opt=it_gdcssomp_opt; end
            %%% Results for channel tracking
            if MODGDCSSOMP_DFT==1, sp.MSE_MODGDCSSOMP_DFT=10*log10(MSE_MODGDCSSOMP_DFT/MSE_all_N); sp.it_modgdcssomp_dft=it_modgdcssomp_dft; end
            if MODGDCSSOMP_OPT==1, sp.MSE_MODGDCSSOMP_OPT=10*log10(MSE_MODGDCSSOMP_OPT/MSE_all_N); sp.it_modgdcssomp_opt=it_modgdcssomp_opt; end
            
            
%             cd results                        %% in case the results are to be stored in a separate directory
%             tmpl_s='results_CCE_%2.4d';       %% file name for saving the results
%             fname=sprintf(tmpl_s,iter);       %% in case you want to save each iteration separately
            fname='results_CCE';
            save(fname, 'sp');
%             cd ..
        end  %% steps
    end  %% reln
end  %% iter