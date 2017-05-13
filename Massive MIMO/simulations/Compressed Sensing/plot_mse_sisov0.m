


figura = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[EbNoVec(1)-0.01, EbNoVec(end)],'ylim',[1e-8 1000]);
xlabel('Eb/No (dB)'); ylabel('MSE'); set(figura,'NumberTitle','off');
set(figura, 'renderer', 'zbuffer'); set(figura,'Name','OFDM modulated with QPSK Massive MU-MIMO System');
strTitle = sprintf('Massive MU-MIMO Channel Estimation on Uplink - Np: %d',Np);
title(strTitle);


teste2(1,:) = mse_siso(70,1,:);

semilogy(EbNoVec,teste2,'k')
semilogy(EbNoVec, est_error_ls, 'bs');
semilogy(EbNoVec, est_error_mmse, 'k*');
semilogy(EbNoVec, est_error_omp, 'ro');
legend('MSE-SISO','LS','MMSE','OMP');
hold off