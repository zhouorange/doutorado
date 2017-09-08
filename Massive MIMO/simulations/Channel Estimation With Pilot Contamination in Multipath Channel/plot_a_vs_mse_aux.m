clear all;close all;clc

load('.mat');
theoretical_mmse_error_m20 = theoretical_mmse_error;
theoretical_ls_error_m20 = theoretical_ls_error;
prop_4_ana_error_vec_m20 = prop_4_ana_error_vec;

clear 'theoretical_mmse_error' 'theoretical_ls_error' 'prop_4_ana_error_vec'

load('.mat');
theoretical_mmse_error_m30 = theoretical_mmse_error;
theoretical_ls_error_m30 = theoretical_ls_error;
prop_4_ana_error_vec_m30 = prop_4_ana_error_vec;

clear 'theoretical_mmse_error' 'theoretical_ls_error' 'prop_4_ana_error_vec'

load('.mat');
theoretical_mmse_error_m90 = theoretical_mmse_error;
theoretical_ls_error_m90 = theoretical_ls_error;
prop_4_ana_error_vec_m90 = prop_4_ana_error_vec;

clear 'theoretical_mmse_error' 'theoretical_ls_error' 'prop_4_ana_error_vec'

% fdee_figure = figure;
% subplot(1,2,1)
% loglog(a,theoretical_ls_error_m30,'--b','MarkerSize',7);
% hold on;
% loglog(a,theoretical_mmse_error_m30,'--r');
% loglog(a,real(prop_4_ana_error_vec_m30),'--c','MarkerSize',7);
% hold off
% grid on;
% xlabel('a')
% ylabel('MSE')
% legend('LS (ana)','MMSE (ideal)','Prop. (ana)');

subplot(1,2,1)
loglog(a,theoretical_ls_error_m90,'--b*','MarkerSize',7);
hold on;
loglog(a,theoretical_mmse_error_m90,'--r*');
loglog(a,real(prop_4_ana_error_vec_m90),'--c*','MarkerSize',7);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('LS (ana)','MMSE (ideal)','Prop. (ana)');

subplot(1,2,2)
loglog(a,theoretical_ls_error_m20,'--b*','MarkerSize',7);
hold on;
loglog(a,theoretical_mmse_error_m20,'--r*');
loglog(a,real(prop_4_ana_error_vec_m20),'--c*','MarkerSize',7);
hold off
grid on;
xlabel('a')
ylabel('MSE')
legend('LS (ana)','MMSE (ideal)','Prop. (ana)');