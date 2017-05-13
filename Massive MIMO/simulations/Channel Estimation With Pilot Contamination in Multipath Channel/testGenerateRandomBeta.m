clear all;close all;clc;

L = 7;                                                              % Number of cells.
K = 10;                                                             % Number of single-antenna terminals in each cell, i.e., number of transmitt antennas.
M = 100;                                                             % Number of antennas at the base station of each cell (In this case, i.e., uplink, it gives the number of receiving antennas).

cellRadius = 1000;                                                  % Radius given in meters.
cellHole = 100;                                                     % Cell hole in meters.

% Large scale fading.
sshadow = 3;                                                        % Shadow-fading standard deviation in dB.
gamma = 2.8;

[distances, angles] = GenerateRandomBetav2(L,K,cellRadius,cellHole,sshadow,gamma);

figure;

for n=0:1:L-1
    
    if(n==0)
        offsetX = 0;
        offsetY = 0;
        plot(offsetX, offsetY, 'rs', 'MarkerFaceColor',[1,0,0], 'MarkerEdgeColor',[1 0 0]); % eNodeB location.
        hold on;
    else
        offsetX = 2*cellRadius*cos((2*pi*(n-1))/(L-1));
        offsetY = 2*cellRadius*sin((2*pi*(n-1))/(L-1));
        plot(offsetX, offsetY, 'rs', 'MarkerFaceColor',[1,0,0], 'MarkerEdgeColor',[1 0 0]); % eNodeB location.
    end
    
    Circlex = offsetX + cellRadius*cos(2*pi*(0:999)/1000);
    Circley = offsetY + cellRadius*sin(2*pi*(0:999)/1000);
    plot(Circlex.', Circley.', 'k');
    
    Circlex = offsetX + cellHole*cos(2*pi*(0:999)/1000);
    Circley = offsetY + cellHole*sin(2*pi*(0:999)/1000);
    plot(Circlex.', Circley.', 'r');
    
    UEx = distances(n+1,:).*cos(angles(n+1,:));
    UEy = distances(n+1,:).*sin(angles(n+1,:));
    plot(UEx, UEy, 'b*');
end

grid on;
hold off;