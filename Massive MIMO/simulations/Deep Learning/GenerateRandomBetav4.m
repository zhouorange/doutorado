function beta_i = GenerateRandomBetav4(L,K,cellRadius,cellHole,sshadow,gamma)

% **************************** Cell Layout. *******************************

distances = zeros(L,K);
angles = zeros(L,K);
beta_i = zeros(L,K);

% Generate the UEs locations for all cells.
for n=0:1:L-1
    
    if(n==0) % Position UEs in the main cell.
        offsetX = 0;
        offsetY = 0;
    else
        offsetX = 2*cellRadius*cos((2*pi*(n-1))/(L-1));
        offsetY = 2*cellRadius*sin((2*pi*(n-1))/(L-1));
    end
    
    % Generate a random radius value within the range cellHole to cellRadius for each one of the terminals.
    radius = cellHole + (cellRadius-cellHole).*rand(1,K);
    
    % Generate an random angle value within the range 0 to 2*pi for each one of the terminals.
    angle = 2*pi*rand(1,K);
    
    UEx = offsetX + radius.*cos(angle);
    UEy = offsetY + radius.*sin(angle);
    
    % Calculate distance to cell i (main cell).
    distances(n+1,:) = sqrt(UEx.^2 + UEy.^2);
    angles(n+1,:) = atan2(UEy,UEx);
    
    psi_db = sqrt(sshadow)*randn(1,K);
    psi_lin = 10.^(psi_db./10);
    
    beta_i(n+1,:) = psi_lin./((distances(n+1,:)./cellHole).^gamma);
end
