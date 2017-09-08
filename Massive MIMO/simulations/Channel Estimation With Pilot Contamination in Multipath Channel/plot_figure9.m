clear all;close all;clc

SNR = -50:4:30;                           % Signal-to-Noise ratio, given in dB.

M = [10 30 100 200 500];                             % Number of antennas.
K = 10;                             % Number of single-antenna users.
P = 20;                             % Channel Length (Finite Impulse Response - FIR).
L = 7;                              % Number of cells.

N = getPilotLength(K,P);            % Pilot length is set according to K and P.
q = 10.^(SNR./10);                  % Uplink pilot power.

a = 0.05;                           % Cross-cell interference level.

isFixedValue = true;
if(isFixedValue)
    beta = a*ones(L,K);
    for ll=1:1:L
        for kk=1:1:K
            if(kk==ll)
                beta(ll,kk) = 1;
            end
        end
    end
else
    beta = abs(randn(L,K));
end
k_idx = 1;
beta_sum = 0;
for l_idx=1:1:L
    beta_sum = beta_sum + beta(l_idx,k_idx);
end
beta111 = beta(1,1);

% Simulation starts here.
distance = zeros(length(M),length(q));
for m_idx=1:1:length(M)
    for q_idx = 1:1:length(q)
        epsilon_11 = (beta_sum + 1/(q(q_idx)*N));
        
        distance(m_idx,q_idx) = (1./((M(m_idx)*P)-1)).*(beta111.^2).*(1./epsilon_11);
    end
    fprintf(1, 'M: %d\n',M(m_idx));
end

fontSize = 12;
lineWidth = 1;
markerSize = 8;

semilogy(SNR,distance(1,:),'s-','MarkerSize', markerSize, 'LineWidth', lineWidth);
hold on
semilogy(SNR,distance(2,:),'o-','MarkerSize', markerSize, 'LineWidth', lineWidth);
semilogy(SNR,distance(3,:),'x-','MarkerSize', markerSize, 'LineWidth', lineWidth);
semilogy(SNR,distance(4,:),'^-','MarkerSize', markerSize, 'LineWidth', lineWidth);
semilogy(SNR,distance(5,:),'+-','MarkerSize', markerSize, 'LineWidth', lineWidth);
xlabel('SNR[dB]')
ylabel('Distance between Prop. and MMSE estimators')
%axis([SNR(1) SNR(length(SNR)) 9e-6 0.005])
leg1 = legend('M = 10','M = 30','M = 100','M = 200','M = 500','Location','southeast');
grid on
hold off

set(gca, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(leg1, 'FontName', 'Times New Roman', 'FontSize', fontSize)
set(gca, 'defaultAxesFontName', 'Times New Roman')
set(gca, 'defaultTextFontName', 'Times New Roman')

scaleFactor = 1.6;
set(gcf, 'Position', [100, 100, ceil(scaleFactor*560), ceil(scaleFactor*420)])

