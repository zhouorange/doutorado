%%% This file sets up all receiver related data
%%%
%%% The following fields have to be set up here:
%%% Rx.G                 Array data (using make_geometry function)
%%% Rx.array_rot[x|y|z]  Array rotation angles
%%% Rx.positionc         Position of the array center
%%%
%%% NOTE: All paramters have to be fields of the struct Geometry
%%%
%%% modified by Daniel Eiwen, 2009-2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state', sum(100*clock));


% for easier access
T = Geometry.T;

load('parameters');     %% load the parameters N,L,steps,NT,NR,bandwidth,f0,LL

% --- set array type and array-related data ---
% Rx.G = make_array(array_type, antenna_number, antenna_spacing)
% (see "Available Arrays" in Help menu)
Rx.G = make_array('ULA', NR, Geometry.lambda0/2,'vertical_dipole');  % ULA with NR antennas

% --- set array rotation for all time snapshots ---
% rotx rotates around the x-axis, roty rotates around the y-axis, 
% rotz rotates around the z-axis
% all angles are in radians
Rx.array_rotx = 0*ones(1,T);
Rx.array_roty = 0*ones(1,T);
Rx.array_rotz = 0*ones(1,T);

% --- set position of the Tx for all time snapshots ---
v_ms=50*rand(1);          % random velocity of up to 50 m/s
v = rand(3,1)-0.5;        % random direction vector [x y z]'
v = v_ms*v/norm(v);    
a_ms=7*rand(1);           % random acceleration of up to 7 m/s^2
a = rand(3,1)-0.5;        % random direction vector [x y z]'
a = a_ms*a/norm(a);    
startpoint = [0 0 2]';    % start point in [x y z]'

Rx.positionc = zeros(3,T);  % allocate memory
Rx.positionc(1,:) = startpoint(1)+v(1)*Geometry.time+a(1)/2*Geometry.time.^2;
Rx.positionc(2,:) = startpoint(2)+v(2)*Geometry.time+a(2)/2*Geometry.time.^2;
Rx.positionc(3,:) = startpoint(3)+v(3)*Geometry.time+a(3)/2*Geometry.time.^2;

% --- assign the Tx structure as first user of the Geometry struct ---
Geometry.Rx = Rx;
