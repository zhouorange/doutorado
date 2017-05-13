%%% This file generates the information about the environment
%%%
%%% It has to define the following fields:
%%% S.positions
%%% S.coefficients
%%% paths
%%% Obstacles.[x|y|z]
%%%
%%% NOTE: All paramters have to be fields of the struct Geometry,
%%%
%%% modified by Daniel Eiwen, 2009-2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- set the scatterer postions ---
% can be defined in an arbitrary way - the function scatterize can be used
% for clusters

rand('state', sum(100*clock));

T=Geometry.T;       %% for easy access

load('parameters'); %% load the parameters N,L,steps,NT,NR,bandwidth,f0,LL

fact=1.2;           %% factor determining how far "outside" the transmitter the clusters are allowed to be

xmin = 0;           %% minimal/maximal x/y/z coordinate for clusters
xmax = 1500*fact;
ymin = -500;
ymax = 500;
zmin = 0;
zmax = 40;

Numb_Cluster=10;        %% Total number of clusters of scatterers
Numb_Scat=10;           %% Number of scatterers per cluster
Numb_Surr=3;            %% Number of clusters surrounding the receiver
Rad_Surr=100;           %% [m], radius of the "surrounding area" of the receiver

hlp=randc(Numb_Surr,1);
hlp=2*Rad_Surr*(rand(Numb_Surr,1)-0.5).*hlp./abs(hlp);

% --- set position of the clusters at the start ---
X = [real(hlp).',rand(1,Numb_Cluster-Numb_Surr)*abs(xmax-xmin)+min([xmin,xmax])];
Y = [imag(hlp).',rand(1,Numb_Cluster-Numb_Surr)*abs(ymax-ymin)+min([ymin,ymax])];
Z = rand(1,Numb_Cluster)*abs(zmax-zmin)+min([zmin,zmax]);

startpoint = [X;Y;Z]; % start point in [x y z]'

v_ms = rand(1,Numb_Cluster)*50;       % random velocity of up to 50 m/s
v = rand(3,Numb_Cluster)-0.5;         % random direction vector [x y z]'
a_ms = rand(1,Numb_Cluster)*7;        % random acceleration of up to 7 m/s^2
a = rand(3,Numb_Cluster)-0.5;         % random direction vector [x y z]'
for kk=1:Numb_Cluster
    v(:,kk) = v_ms(kk)*v(:,kk)/norm(v(:,kk));
    a(:,kk) = a_ms(kk)*a(:,kk)/norm(a(:,kk));
end

Posi=zeros(3,T,Numb_Cluster);
for kk=1:Numb_Cluster
    Posi(1,:,kk) = startpoint(1,kk)+v(1,kk)*Geometry.time+a(1,kk)/2*Geometry.time.^2;
    Posi(2,:,kk) = startpoint(2,kk)+v(2,kk)*Geometry.time+a(2,kk)/2*Geometry.time.^2;
    Posi(3,:,kk) = startpoint(3,kk)+v(3,kk)*Geometry.time+a(3,kk)/2*Geometry.time.^2;
end


%%% create clusters of Numb_Scat scatterer each
%%% distributed in a sphere of some radius arounf the center points 

%%% clusters 1-3 surround the receiver
% cluster 1
S1.positions = Posi(:,:,1);
S1 = scatterize(S1,'spherical',Numb_Scat,5);
% cluster 2
S2.positions = Posi(:,:,2);
S2 = scatterize(S2,'spherical',Numb_Scat,10);
% cluster 3
S3.positions = Posi(:,:,3);
S3 = scatterize(S3,'spherical',Numb_Scat,20);

% cluster 4
S4.positions = Posi(:,:,4);
S4 = scatterize(S4,'spherical',Numb_Scat,5);
% cluster 5
S5.positions = Posi(:,:,5);
S5 = scatterize(S5,'spherical',Numb_Scat,10);
% cluster 6
S6.positions = Posi(:,:,6);
S6 = scatterize(S6,'spherical',Numb_Scat,20);
% cluster 7
S7.positions = Posi(:,:,7);
S7 = scatterize(S7,'spherical',Numb_Scat,30);
% cluster 8
S8.positions = Posi(:,:,8);
S8 = scatterize(S8,'spherical',Numb_Scat,10);
% cluster 9
S9.positions = Posi(:,:,9);
S9 = scatterize(S9,'spherical',Numb_Scat,20);
% cluster 10
S10.positions = Posi(:,:,10);
S10 = scatterize(S10,'spherical',Numb_Scat,30);


% concatenation
S.positions = cat(3,S1.positions,S2.positions,S3.positions,S4.positions,S5.positions,S6.positions,S7.positions,S8.positions,S9.positions,S10.positions);

% --- set number of scatterers
S.O = size(S.positions,3);

% --- set the scatterer coefficients ---
S.coefficients = repmat(0.3*rand(S.O,1).*exp(1j*rand(S.O,1)*2*pi),[1,T]);

Geometry.S = S;

% --- define the paths matrix to contain all single reflections ---
% NOTE: path must be one of the classes uint8 or uint16
Geometry.paths = paths_addsingles(uint8([]), Geometry.S);

% --- define the Obstacles ---
% the fields define the bounds of the obstacles in the rows (for each
% dimension) and the number of columns is the number of obstacles
% Geometry.Obstacles.x = [24 26].';
% Geometry.Obstacles.y = [-5 -7].';
% Geometry.Obstacles.z = [0   4].';
