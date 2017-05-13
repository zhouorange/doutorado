%%% gen_channel_cce.m
%%%
%%% Generates a MIMO channel adapted to channel tracking using the free toolbox IlmProp, 
%%% which has to be installed and adapted first. 
%%% 
%%% Input   -     N     = length of one OFDM symbol
%%%               L     = number of OFDM symbols
%%%               steps = number of tracking steps
%%%               NT,NR = number of transmit/receive antennas
%%%               bandwidth = bandwidth
%%%               f0    = carrier frequency
%%%               LL    = number of frequency snapshots
%%%               iter  = iteration count, just needed if the geometry and 
%%%                       the channel are to be stored
%%%
%%% Output  -     H     = channel array of size NR x NT x N*L*steps x LL
%%%
%%% Usage   -     H     = gen_channel_cce(N,L,steps,NT,NR,bandwidth,f0,LL,iter)
%%%
%%% COPYRIGHT : (c) Daniel Eiwen, 2012

function H=gen_channel_cce(N,L,steps,NT,NR,bandwidth,f0,LL,iter)
clear H;
clear Geometry;
rand('state', sum(100*clock));

global Geometry
cd IlmProp          %% change to wherever your IlmProp is located
save('parameters','N','L','steps','NT','NR','bandwidth','f0','LL');
Geometry = init_geometry('CCE', 'geometry_scripts');
H = forward_freq(Geometry);
cd ..

%%% save Geometry and channel for later use
if nargin==8
    save(sprintf('geometry_%2.3d',iter),'Geometry');
    save(sprintf('channel_%2.3d',iter),'H');
end