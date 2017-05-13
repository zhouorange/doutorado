 function [distances] = GenerateRandomBeta(L,K,cellRadius,cellHole,sshadow,gamma)

% L = 7;                                                              % Number of cells.
% K = 10;                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
% M = 100;                                                             % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).
% 
% cellRadius = 1000;                                                  % Radius given in meters.
% cellHole = 100;                                                     % Cell hole in meters.
% 
% % Large scale fading.
% sshadow = 8;                                                        % Shadow-fading standard deviation in dB.
% gamma = 3.8;                                                        % Decay exponent: Urban area cellular radio 2.7 - 3.5

% **************************** Cell Layout. *******************************

distances = zeros(L,K);

% ********** Center (Main) Cell  **********

Circlex = cellRadius*cos(2*pi*(0:99)/100);
Circley = cellRadius*sin(2*pi*(0:99)/100);

Circlex = cellHole*cos(2*pi*(0:99)/100);
Circley = cellHole*sin(2*pi*(0:99)/100);

% Position UEs in the main cell.
% Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
radius = cellHole + (cellRadius-cellHole).*rand(1,K);

% Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
angle = 2*pi*rand(1,K);

UEx = radius.*cos(angle);
UEy = radius.*sin(angle);

distances(1,:) = sqrt(UEx.^2 + UEy.^2); % Distance in meters.

% Generate the UEs locations of cells around the main cell.
for n=0:1:L-2
    offsetX = 2*cellRadius*cos((2*pi*n)/(L-1));
    offsetY = 2*cellRadius*sin((2*pi*n)/(L-1));

    Circlex = offsetX + cellRadius*cos(2*pi*(0:99)/100);
    Circley = offsetY + cellRadius*sin(2*pi*(0:99)/100);
    
    Circlex = offsetX + cellHole*cos(2*pi*(0:99)/100);
    Circley = offsetY + cellHole*sin(2*pi*(0:99)/100);
    
    % Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
    radius = cellHole + (cellRadius-cellHole).*rand(1,K);
    
    % Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
    angle = 2*pi*rand(1,K);
    
    UEx = offsetX + radius.*cos(angle);
    UEy = offsetY + radius.*sin(angle);
    
    % Calculate distance to cell i (main cell).
    distances(n+2,:) = sqrt(UEx.^2 + UEy.^2);
end
