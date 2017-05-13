clear all;close all;clc

v = linspace(0.01,340.29,10000);  % Speed in m/s.  
c = 3e8;                        % Light speed in m/s.

fc1 = 1800e6;                   % Carrier Freq. in MHz.
fd_fc1 = (v*fc1)/c;
Tc_fc1 = sqrt((9/16*pi)).*(1./fd_fc1);

fc2 = 800e6;                   % Carrier Freq. in MHz.
fd_fc2 = (v*fc2)/c;
Tc_fc2 = sqrt((9/16*pi)).*(1./fd_fc2);

fc3 = 2100e6;                   % Carrier Freq. in MHz.
fd_fc3 = (v*fc3)/c;
Tc_fc3 = sqrt((9/16*pi)).*(1./fd_fc3);

fc4 = 400e6;                   % Carrier Freq. in MHz.
fd_fc4 = (v*fc4)/c;
Tc_fc4 = sqrt((9/16*pi)).*(1./fd_fc4); 

fc5 = 2600e6;                   % Carrier Freq. in MHz.
fd_fc5 = (v*fc5)/c;
Tc_fc5 = sqrt((9/16*pi)).*(1./fd_fc5);

semilogy(v,Tc_fc4,'MarkerSize',7,'LineWidth',1)
hold on
semilogy(v,Tc_fc2,'MarkerSize',7,'LineWidth',1)
semilogy(v,Tc_fc1,'MarkerSize',7,'LineWidth',1)
semilogy(v,Tc_fc3,'MarkerSize',7,'LineWidth',1)
semilogy(v,Tc_fc5,'MarkerSize',7,'LineWidth',1)
grid on
xlabel('v [m/s]');
ylabel('Tc [s]');
legend('f_{c} = 400 MHz','f_{c} = 800 MHz','f_{c} = 1800 MHz','f_{c} = 2100 MHz','f_{c} = 2600 MHz');
axis([v(1) v(length(v)) 0.0003 50]);
hold off