clear all;close all;clc

plot_fig = true;                                                   % Disable/enable figure plotting.

L = 7;                                                              % Number of cells.
K = 10;                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 100;                                                             % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).

cellRadius = 1000;                                                  % Radius given in meters.
cellHole = 100;                                                     % Cell hole in meters.

% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

% **************************** Cell Layout. *******************************
% Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
radius = cellHole + (cellRadius-cellHole).*rand(1,K);

% Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
angle = 2*pi*rand(1,K);

% Plot the position of each terminal inside the cell.
figure;

% ********** Center (Main) Cell  **********
% Plot cell center, where eNodeB is located.
plot(0, 0, 'rs', 'MarkerFaceColor',[1,0,0], 'MarkerEdgeColor',[1 0 0]);
hold on

Circlex = cellRadius*cos(2*pi*(0:99)/100);
Circley = cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');

Circlex = cellHole*cos(2*pi*(0:99)/100);
Circley = cellHole*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'r');

% ********** 2nd Cell  **********
n = 0;
offsetX = 2*cellRadius*cos((2*pi*n)/(L-1));
offsetY = 2*cellRadius*sin((2*pi*n)/(L-1));
plot(offsetX, offsetY, 'rs', 'MarkerFaceColor',[1,0,0], 'MarkerEdgeColor',[1 0 0]); % eNodeB location.
Circlex = offsetX + cellRadius*cos(2*pi*(0:99)/100);
Circley = offsety + cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');

% ********** 3rd Cell  ********** 
n = 1;
Circlex = 2*cellRadius*cos((2*pi*n)/(L-1)) + cellRadius*cos(2*pi*(0:99)/100);
Circley = 2*cellRadius*sin((2*pi*n)/(L-1)) + cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');

% ********** 4th Cell  ********** 
n = 2;
Circlex = 2*cellRadius*cos((2*pi*n)/(L-1)) + cellRadius*cos(2*pi*(0:99)/100);
Circley = 2*cellRadius*sin((2*pi*n)/(L-1)) + cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');

% ********** 5th Cell  ********** 
n = 3;
Circlex = 2*cellRadius*cos((2*pi*n)/(L-1)) + cellRadius*cos(2*pi*(0:99)/100);
Circley = 2*cellRadius*sin((2*pi*n)/(L-1)) + cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');

% ********** 6th Cell  ********** 
n = 4;
Circlex = 2*cellRadius*cos((2*pi*n)/(L-1)) + cellRadius*cos(2*pi*(0:99)/100);
Circley = 2*cellRadius*sin((2*pi*n)/(L-1)) + cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');

% ********** 7th Cell  ********** 
n = 5;
Circlex = 2*cellRadius*cos((2*pi*n)/(L-1)) + cellRadius*cos(2*pi*(0:99)/100);
Circley = 2*cellRadius*sin((2*pi*n)/(L-1)) + cellRadius*sin(2*pi*(0:99)/100);
plot(Circlex.', Circley.', 'k');



UEx = radius.*cos(angle);
UEy = radius.*sin(angle);
plot(UEx, UEy, 'b*');

grid on;
hold off;

% Calculate path-loss for all users, in meters.
path_loss = radius.^gamma;

% Calculate shadowing for each one of different channels.
shadowing_att = lognrnd(0,sshadow,1,K);
% Large scale fading calculated according to Marzetta.
largeScaleFading = repmat(sqrt(shadowing_att./path_loss), M, 1);
