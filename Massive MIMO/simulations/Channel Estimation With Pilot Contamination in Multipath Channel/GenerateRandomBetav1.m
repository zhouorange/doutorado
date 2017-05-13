 function [distances] = GenerateRandomBetav1(L,K,cellRadius,cellHole,sshadow,gamma)

% **************************** Cell Layout. *******************************

distances = zeros(L,K);

% ********** Center (Main) Cell  **********

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
    
    % Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
    radius = cellHole + (cellRadius-cellHole).*rand(1,K);
    
    % Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
    angle = 2*pi*rand(1,K);
    
    UEx = offsetX + radius.*cos(angle);
    UEy = offsetY + radius.*sin(angle);
    
    % Calculate distance to cell i (main cell).
    distances(n+2,:) = sqrt(UEx.^2 + UEy.^2);
end
