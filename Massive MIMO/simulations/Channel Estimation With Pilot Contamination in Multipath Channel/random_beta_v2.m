clear all;close all;clc

L = 7;                                                              % Number of cells.
K = 10;                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 100;                                                             % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).

cellRadius = 1000;                                                  % Radius given in meters.
cellHole = 100;                                                     % Cell hole in meters.

% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

% **************************** Cell Layout. *******************************

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

% Position UEs in the main cell.
% Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
radius = cellHole + (cellRadius-cellHole).*rand(1,K);

% Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
angle = 2*pi*rand(1,K);

UEx = radius.*cos(angle);
UEy = radius.*sin(angle);
plot(UEx, UEy, 'b*');

for n=0:1:L-2
    offsetX = 2*cellRadius*cos((2*pi*n)/(L-1));
    offsetY = 2*cellRadius*sin((2*pi*n)/(L-1));
    plot(offsetX, offsetY, 'rs', 'MarkerFaceColor',[1,0,0], 'MarkerEdgeColor',[1 0 0]); % eNodeB location.
    Circlex = offsetX + cellRadius*cos(2*pi*(0:99)/100);
    Circley = offsetY + cellRadius*sin(2*pi*(0:99)/100);
    plot(Circlex.', Circley.', 'k');
    
    Circlex = offsetX + cellHole*cos(2*pi*(0:99)/100);
    Circley = offsetY + cellHole*sin(2*pi*(0:99)/100);
    plot(Circlex.', Circley.', 'r');
    
    % Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
    radius = cellHole + (cellRadius-cellHole).*rand(1,K);
    
    % Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
    angle = 2*pi*rand(1,K);
    
    UEx = offsetX + radius.*cos(angle);
    UEy = offsetY + radius.*sin(angle);
    plot(UEx, UEy, 'b*');
    
    % Calculate distance to cell i (main cell).
    distances = sqrt(UEx.^2 + UEy.^2);
    angles = atan(UEy./UEx);
end

grid on;
hold off;

% Calculate path-loss for all users, in meters.
path_loss = radius.^gamma;

% Calculate shadowing for each one of different channels.
shadowing_att = lognrnd(0,sshadow,1,K);
% Large scale fading calculated according to Marzetta.
largeScaleFading = repmat(sqrt(shadowing_att./path_loss), M, 1);
