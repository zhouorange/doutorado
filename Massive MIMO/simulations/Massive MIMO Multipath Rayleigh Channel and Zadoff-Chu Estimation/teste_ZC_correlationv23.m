clear all;close all;clc

% ------- Definitions -------
u = 129;
Nzc = 839;
NIFFT = 2048;
NIDFT = 24576;
NSUBFRAME = 30720;
Ncp = 3168;
NG = 2976;
v = [0];
Ncs = 13;

% ------- Generate Root Zadoff-Chu sequence. -------
n = [0:1:(Nzc-1)];
xu_root = exp(-1i*(pi*u.*n.*(n+1))./Nzc);

% ------- Apply Cyclic Shift to Root Zadoff-Chu sequence. -------
Cv = v(1)*Ncs;
xuv = xu_root(mod((n+Cv),Nzc)+1);

zu = complex(zeros(1,Nzc),zeros(1,Nzc));
for l=0:1:Nzc-1
    for n=0:1:Nzc-1
        zu(l+1) = zu(l+1) + xuv(n+1)*conj(xu_root(mod((n+l),Nzc)+1));
    end
end

% Noise
w = wgn(1,Nzc,0,'complex');
zu = complex(zeros(1,Nzc),zeros(1,Nzc));
for l=0:1:Nzc-1
    for n=0:1:Nzc-1
        zu(l+1) = zu(l+1) + w(n+1)*conj(xu_root(mod((n+l),Nzc)+1));
    end
end
noise = abs(zu/Nzc).^2;

zu = complex(zeros(1,Nzc),zeros(1,Nzc));
y = xuv+w;
for l=0:1:Nzc-1
    for n=0:1:Nzc-1
        zu(l+1) = zu(l+1) + y(n+1)*conj(xu_root(mod((n+l),Nzc)+1));
    end
end
zu = zu/Nzc;

aux = zu(1) - complex(1,0);

zu(1) = zu(1) - aux;




a=1;