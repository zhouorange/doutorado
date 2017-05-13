clear all;clc


% We start by defining some common simulation parameters
N = 1000000;                  % Number of transmit antennas
modOrd = 2;             % constellation size = 2^modOrd

%%
% and set up the simulation.

% Create a local random stream to be used by random number generators for
% repeatability.
hStr = RandStream('mt19937ar');

% Create PSK modulator and demodulator System objects
hMod   = comm.PSKModulator(...
    'ModulationOrder',  2^modOrd, ...
    'PhaseOffset',      0, ...
    'BitInput',         true);

% Create random bit vector to modulate
msg = randi(hStr, [0 1], [N*modOrd, 1]);

% Modulate data
txSig = step(hMod, msg);

aaa=var(txSig.');

a=1;