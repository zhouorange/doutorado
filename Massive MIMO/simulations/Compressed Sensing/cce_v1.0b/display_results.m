%%% display_results.M
%%%
%%% displays the results of CompChaEst.m
%%% 
%%% COPYRIGHT : (c) Daniel Eiwen, 2012
%%%
%%% Note that here we just show HOW TO display the results. Adapt this file 
%%% to display the specific results you are interested in. 
%%% Note also that in the case of group size (1,1) all the algorithms
%%% taking group-sparsity into account coincide with their conventional
%%% counterparts, so to compare them compare the group sizes different from
%%% (1,1) to group size (1,1). Furthermore, MOD-GDCSSOMP coincides with
%%% conventional GDCSSOMP in the case of suppsize=0, i.e. if no support
%%% information is given.

load('results_CCE');

ginum=length(sp.gr_i);
gmnum=length(sp.gr_m);

MSE_SISO=cell(ginum,gmnum);
MSE_MIMO=cell(ginum,gmnum);
MSE_Track_DFT=cell(ginum,gmnum);
MSE_Track_OPT=cell(ginum,gmnum);

for gi=1:ginum
    for gm=1:gmnum
        %%% results for channel-per-channel estimation
        if sp.LS==1
            MSE_SISO{gi,gm}=[MSE_SISO{gi,gm},squeeze(sp.MSE_LS)]; 
        end
        if sp.GBPDN_DFT==1
            MSE_SISO{gi,gm}=[MSE_SISO{gi,gm},squeeze(sp.MSE_GBPDN_DFT(gi,gm,:))]; 
        end
        if sp.GCOSAMP_DFT==1
            MSE_SISO{gi,gm}=[MSE_SISO{gi,gm},squeeze(sp.MSE_GCOSAMP_DFT(gi,gm,:))]; 
        end
        if sp.GOMP_DFT==1
            MSE_SISO{gi,gm}=[MSE_SISO{gi,gm},squeeze(sp.MSE_GOMP_DFT(gi,gm,:))]; 
        end
        if sp.GBPDN_OPT==1
            MSE_SISO{gi,gm}=[MSE_SISO{gi,gm},squeeze(sp.MSE_GBPDN_OPT(gi,gm,:))]; 
        end
        if sp.GCOSAMP_OPT==1
            MSE_SISO{gi,gm}=[MSE_SISO{gi,gm},squeeze(sp.MSE_GCOSAMP_OPT(gi,gm,:))]; 
        end
        if sp.GOMP_OPT==1
            MSE_SISO{gi,gm}=[MSE_SISO{gi,gm},squeeze(sp.MSE_GOMP_OPT(gi,gm,:))]; 
        end

        %%% results for multichannel estimation
        if sp.MGBPDN_DFT==1
            MSE_MIMO{gi,gm}=[MSE_MIMO{gi,gm},squeeze(sp.MSE_MGBPDN_DFT(gi,gm,:))]; 
        end
        if sp.MGCOSAMP_DFT==1
            MSE_MIMO{gi,gm}=[MSE_MIMO{gi,gm},squeeze(sp.MSE_MGCOSAMP_DFT(gi,gm,:))]; 
        end
        if sp.GDCSSOMP_DFT==1
            MSE_MIMO{gi,gm}=[MSE_MIMO{gi,gm},squeeze(sp.MSE_GDCSSOMP_DFT(gi,gm,:))]; 
        end
        if sp.MGBPDN_OPT==1
            MSE_MIMO{gi,gm}=[MSE_MIMO{gi,gm},squeeze(sp.MSE_MGBPDN_OPT(gi,gm,:))]; 
        end
        if sp.MGCOSAMP_OPT==1
            MSE_MIMO{gi,gm}=[MSE_MIMO{gi,gm},squeeze(sp.MSE_MGCOSAMP_OPT(gi,gm,:))]; 
        end
        if sp.GDCSSOMP_OPT==1
            MSE_MIMO{gi,gm}=[MSE_MIMO{gi,gm},squeeze(sp.MSE_GDCSSOMP_OPT(gi,gm,:))]; 
        end
        
        %%% results for channel tracking
        if sp.MODGDCSSOMP_DFT==1
            MSE_Track_DFT{gi,gm}=[MSE_Track_DFT{gi,gm},squeeze(sp.MSE_MODGDCSSOMP_DFT(:,gi,gm,:))]; 
        end
        if sp.MODGDCSSOMP_OPT==1
            MSE_Track_OPT{gi,gm}=[MSE_Track_OPT{gi,gm},squeeze(sp.MSE_MODGDCSSOMP_OPT(:,gi,gm,:))]; 
        end
    end
end

%%% legends for the plots
leg_SISO={'LS';'G-BPDN DFT';'G-COSAMP DFT';'G-OMP DFT';'G-BPDN OPT'; 'G-COSAMP OPT';'G-OMP OPT'};
%%% set the entries of alg_SISO corresponding to the algorithms, the results
%%% of which you want to plot, to 1, and the others to 0.
alg_SISO=[sp.LS;sp.GBPDN_DFT;sp.GCOSAMP_DFT;sp.GOMP_DFT;sp.GBPDN_OPT;sp.GCOSAMP_OPT;sp.GOMP_OPT];

leg_MIMO={'Multichannel G-BPDN DFT';'Multichannel G-COSAMP DFT';'G-DCSSOMP DFT';'Multichannel G-BPDN OPT'; 'Multichannel G-COSAMP OPT';'G-DCSSOMP OPT'};
alg_MIMO=[sp.MGBPDN_DFT;sp.MGCOSAMP_DFT;sp.GDCSSOMP_DFT;sp.MGBPDN_OPT;sp.MGCOSAMP_OPT;sp.GDCSSOMP_OPT];

leg_Track={'Modified G-DCSSOMP DFT';'Modified G-DCSSOMP OPT'};
alg_Track=[sp.MODGDCSSOMP_DFT;sp.MODGDCSSOMP_OPT];

%%% Compare the results from channel estimation for the different
%%% channel-per-channel and multichannel methods for each group size
%%% individually
for gi=1:ginum
    for gm=1:gmnum
        figure;
        plot(-10*log10(sp.NOISE_LEVEL),[MSE_SISO{gi,gm},MSE_MIMO{gi,gm}]);
        xlabel('SNR [dB]');
        ylabel('MSE [dB]');
        legend(leg_SISO{alg_SISO==1},leg_MIMO{alg_MIMO==1});
        grid;
        title(sprintf('MSE with group size (%2.1d,%2.1d)',2^(sp.gr_i(gi)),2^(sp.gr_m(gm))));
    end
end

%%% Compare the results from channel estimation for the different group
%%% sizes using only one sparse recovery method (here for example GDCSSOMP
%%% using the DFT basis)
figure;
plot(-10*log10(sp.NOISE_LEVEL),reshape(sp.MSE_GDCSSOMP_DFT,ginum*gmnum,length(sp.NOISE_LEVEL)));
xlabel('SNR [dB]');
ylabel('MSE [dB]');
leggr=cell(ginum+gmnum,1);
for gi=1:ginum
    for gm=1:gmnum
        leggr{(gm-1)*ginum+gi}=sprintf('(%2.1d,%2.1d)',2^(sp.gr_i(gi)),2^(sp.gr_m(gm)));
    end
end
grid;
legend(leggr);
title('MSE for GDCSSOMP for different group sizes');

%%% Show the results for channel tracking for the different support sizes
%%% for each group size individually
if sp.MODGDCSSOMP_DFT+sp.MODGDCSSOMP_OPT>0
    for gi=1:ginum
        for gm=1:gmnum
            figure;
            plot(-10*log10(sp.NOISE_LEVEL),[MSE_Track_DFT{gi,gm}.',MSE_Track_OPT{gi,gm}.']);
            xlabel('SNR [dB]');
            ylabel('MSE [dB]');
            legtr=cell(length(sp.suppsize),1);
            for susi=1:length(sp.suppsize)
                legtr{susi}=sprintf('%2.2d %%',sp.suppsize(susi)*100);
            end
            grid;
            legend(legtr);
            title(sprintf('MSE for channel tracking using MODDCSSOMP for group size (%2.1d,%2.1d) \n for different sizes of the known support',2^(sp.gr_i(gi)),2^(sp.gr_m(gm))));
        end
    end
end

%%% Show the results for channel tracking using one algorithm versus the 
%%% size of the support for the different values of the SNR (again for each
%%% group size individually)
if sp.MODGDCSSOMP_DFT+sp.MODGDCSSOMP_OPT>0
    for gi=1:ginum
        for gm=1:gmnum
            figure;
            plot(100*sp.suppsize,[MSE_Track_DFT{gi,gm},MSE_Track_OPT{gi,gm}]);
            xlabel('Size of the known support [%]');
            ylabel('MSE [dB]');
            legno=cell(length(sp.NOISE_LEVEL),1);
            for reln=1:length(sp.NOISE_LEVEL)
                legno{reln}=sprintf('SNR = %2.2d dB',round(-10*log10(sp.NOISE_LEVEL(reln))));
            end
            grid;
            legend(legno);
            title(sprintf('MSE for channel tracking using MODDCSSOMP for group size (%2.1d,%2.1d) \n versus the size of the known support',2^(sp.gr_i(gi)),2^(sp.gr_m(gm))));
        end
    end
end