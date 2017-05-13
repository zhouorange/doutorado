clear all;close all;clc

v = linspace(0.01,7.2,10000);  % Speed in m/s.  
c = 3e8;                        % Light speed in m/s.

fc1 = 28e9;                   % Carrier Freq. in GHz.
fd_fc1 = (v*fc1)/c;
Tc_fc1 = sqrt((9/16*pi)).*(1./fd_fc1);

fc2 = 38e9;                   % Carrier Freq. in GHz.
fd_fc2 = (v*fc2)/c;
Tc_fc2 = sqrt((9/16*pi)).*(1./fd_fc2);

fc3 = 60e9;                   % Carrier Freq. in GHz.
fd_fc3 = (v*fc3)/c;
Tc_fc3 = sqrt((9/16*pi)).*(1./fd_fc3);

fc4 = 73e9;                   % Carrier Freq. in GHz.
fd_fc4 = (v*fc4)/c;
Tc_fc4 = sqrt((9/16*pi)).*(1./fd_fc4);

semilogy(v,Tc_fc1,'MarkerSize',7,'LineWidth',1)
hold on
semilogy(v,Tc_fc2,'MarkerSize',7,'LineWidth',1)
semilogy(v,Tc_fc3,'MarkerSize',7,'LineWidth',1)
semilogy(v,Tc_fc4,'MarkerSize',7,'LineWidth',1)
semilogy(v,2e-3*ones(1,length(v)),'k--','MarkerSize',7,'LineWidth',1);
grid on
xlabel('v [m/s]');
ylabel('Tc [s]');
legend('f_{c} = 28 GHz','f_{c} = 38 GHz','f_{c} = 60 GHz','f_{c} = 73 GHz','Slot length = 2 ms');
axis([v(1) v(length(v)) Tc_fc1(length(Tc_fc1))-Tc_fc1(length(Tc_fc1))*(70/100) Tc_fc1(1)]);
hold off