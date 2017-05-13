%%% This file sets up all transmitter related date (also multiple user)
%%% 
%%% The following parameters must be set up, where u is the number of users
%%% Tx(1:u).G                 Array data (using make_geometry function)
%%% Tx(1:u).array_rot[x|y|z]  Array rotation angles
%%% Tx(1:u).positionc         Position of the array center
%%% los_coefficients          Coefficients for the LOSs
%%%
%%% NOTE: All paramters have to be fields of the struct Geometry,
%%% for multi-user the Tx structs have to set up completely before
%%% adding them to the struct array Geometry.Tx(:)
%%%
%%% modified by Daniel Eiwen, 2009-2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state', sum(100*clock));

% for easier access
T = Geometry.T;

load('parameters');     %% load the parameters N,L,steps,NT,NR,bandwidth,f0,LL

%%%%%%%%%%%%%%%%%%%%%%%% USER 1 %%%%%%%%%%%%%%%%%%%%%%%

% --- set array type and array-related data ---
% Tx.G = make_array(array_type, antenna_number, antenna_spacing)
% (see "Available Arrays" in Help menu)
Tx.G = make_array ('UCA', NT, Geometry.lambda0/2,'vertical_dipole');  % UCA with NR antennas

% --- set array rotation for all time snapshots ---
% rotx rotates around the x-axis, roty rotates around the y-axis, 
% rotz rotates around the z-axis
% all angles are in radians
Tx.array_rotx = 2*pi*rand*ones(1,T);
Tx.array_roty = 2*pi*rand*ones(1,T);
Tx.array_rotz = 2*pi*rand*ones(1,T);


% --- set position of the Rx for all time snapshots ---
pos = [1500 0 20]';   % pos in x y z [m]
Tx.positionc = repmat(pos, [1,T]);    % fixed Rx for all time snapshots

Geometry.Tx=Tx;

%%%%%%%%%%%%%%%%%%%%%%%% USER 2 %%%%%%%%%%%%%%%%%%%%%%%
% to add a second user, copy paste the first, change the parameters, 
% and then assign:
% Geometry.Tx(2) = Tx;

%%%%%%%%%%%%%%%%%%%%% LOS COEFFS %%%%%%%%%%%%%%%

% --- set los coefficients ---
Geometry.los_coefficients = 0*ones(length(Geometry.Tx), T);
% artificially setting the LOS to 0. This could be done also by adding an
% obstacle properly positioned. In this way it's easier though. 