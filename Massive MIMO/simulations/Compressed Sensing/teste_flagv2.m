clear all;close all;clc

NFFT = 2048;

Np = 16;
K = 10;

NCP = 512;

source = randi(NFFT,K,NFFT-K*Np);

%% @@@@@@@@@@@@@@@@ Create Pilots @@@@@@@@@@@@@@@@@@@@@@@@@@@@
fo = 1000; % in Hz.                             % Pilot frequency.
P = exp(1i*2*pi*fo*(0:1:Np-1)/Np);              % Normalized pilot signal, i.e., unit power.

Pdiag = diag(P);

F = fft(eye(NFFT));
F = F(:,1:NCP);


%------------ Create Pilots ------------
flag_data = false(K,NFFT);
flag_pilots = false(K,NFFT);

ppos = 0;
while(length(unique(ppos)) ~= Np*K)
    ppos = randi(NFFT,1,Np*K);
end
ppos = reshape(ppos,K,Np);
ppos = sort(ppos,2);

for l_idx=1:1:K
    for c_idx=1:1:Np
        flag_data(:,ppos(l_idx,c_idx)) = true(K,1);
        flag_pilots(l_idx,ppos(l_idx,c_idx)) = true;
    end
end

% Split source among K terminals.
Tx = complex(zeros(K,NFFT),zeros(K,NFFT));
Tx(~flag_data) = reshape(source, K, numel(source)/K);
Tx(flag_pilots) = repmat(P,K,1);
