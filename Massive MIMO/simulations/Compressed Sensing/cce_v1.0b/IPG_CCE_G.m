%%% This file sets up the general geometry paramters
%%%
%%% It gets called first in geometry generation. The following paramters have
%%% to be set in this file:
%%% name, notes, T, time, F, freq and either f0 or lambda0 
%%%
%%% NOTE: All paramters have to be fields of the struct Geometry
%%%
%%% modified by Daniel Eiwen, 2009-2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- set Geometry name and some notes ---
Geometry.name = 'tracking';
Geometry.notes = 'Sparse Channel';

load('parameters');     % load the parameters N,L,steps,NT,NR,bandwidth,f0,LL

% --- set number of time snapshots ---
Geometry.T = steps*N*L;

% --- set time axis (in seconds) ---
end_time = steps*N*L/bandwidth; % [sec]
Geometry.time = linspace(0, end_time, Geometry.T);

% --- set either center frequency (Hz) or center wavelength (m)
Geometry.f0 = f0; 

% --- set number of frequency bins ---
Geometry.F = LL;

% --- set frequency axis in pass band (Hz) ---
fscan = linspace(-bandwidth/2, bandwidth/2, Geometry.F);
Geometry.freq = Geometry.f0+fscan;