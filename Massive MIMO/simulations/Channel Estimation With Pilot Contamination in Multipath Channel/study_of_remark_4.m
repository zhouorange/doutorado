clear all;close all;clc

SNR = 100;                           % Signal-to-Noise ratio, given in dB.

M = 30; %10:20:100;                             % Number of antennas.
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

semilogy(SNR,distance(1,:));
hold on
semilogy(SNR,distance(2,:));
semilogy(SNR,distance(3,:));
semilogy(SNR,distance(4,:));
semilogy(SNR,distance(5,:));
xlabel('SNR[dB]')
ylabel('Distance between Prop. and MMSE estimators')
axis([SNR(1) SNR(length(SNR)) 9e-6 0.005])
legend('M = 10','M = 30','M = 50','M = 70','M = 90','Location','southeast')
grid on
hold off
