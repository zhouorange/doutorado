
load('channel_estimation_mse_vs_num_antennas_K10_P20_L7_N223_v2_20160802T204318.mat')

% Set up a figure for visualizing BER results.
fdee_figure = figure; grid on; hold on;
set(gca,'yscale','log','xlim',[M(1) M(end)],'ylim',[0.23 0.3615]);
xlabel('M'); ylabel('MSE'); set(fdee_figure,'NumberTitle','off');
set(fdee_figure, 'renderer', 'zbuffer');

% Plot results.
semilogy(M(1:m_idx),real(mmse_error_vec(1:m_idx)),'r*', ...
    M(1:m_idx),real(ls_error_vec(1:m_idx)),'b*', ...
    M(1:m_idx),real(prop_1_error_vec(1:m_idx)),'k*', ...
    M(1:m_idx),real(prop_2_error_vec(1:m_idx)),'kv', ...
    M(1:m_idx),real(prop_3_error_vec(1:m_idx)),'ko', ...
    M(1:m_idx),real(prop_4_error_vec(1:m_idx)),'ks' ,'MarkerSize',7)
legend('MMSE (ideal)', 'LS (ideal)', 'Prop. avg=1', 'Prop. avg=5', 'Prop. avg=10', 'Prop. avg=20');

semilogy(M,real(mmse_error_vec),'r-', ...
    M,real(ls_error_vec),'b-', ...
    M,real(prop_1_error_vec),'k-', ...
    M,real(prop_2_error_vec),'k-', ...
    M,real(prop_3_error_vec),'k-', ...
    M,real(prop_4_error_vec),'k-' ,'MarkerSize',7)
hold off;

% % Get timestamp for saving files.
% timeStamp = datestr(now,30);
% fileName = sprintf('channel_estimation_mse_vs_num_antennas_K%d_P%d_L%d_N%d_v0_%s.fig',K,P,L,N,timeStamp);
% savefig(fdee_figure,fileName);
%
% % Save workspace to .mat file.
% clear fdee_figure
% fileName = sprintf('channel_estimation_mse_vs_num_antennas_K%d_P%d_L%d_N%d_v0_%s.mat',K,P,L,N,timeStamp);
% save(fileName);