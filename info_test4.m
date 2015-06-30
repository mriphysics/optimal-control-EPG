% input data for test 4
%% ================ INITIALIZE PARAMETERS =================================
ESP = 2.9 ;%<- echo spacing, ms
T1 = 500; %<- T1,ms
T2 = 400;  %<- T2,ms
klim = 25; % order of maximum k, use 25
max_a = 1.5; % maximum flip = max_a*pi;
red_b1 = 1; % spatial undersampling factor for B1 maps
TSE_factor = 113; % actual length of sequence;
a0 = pi/180*[90,150.66, 90.18,67.32,59.94*ones(1,48)]; % target sequence 
% define piecewise constant indices
singles = 22; % number of individual pulses at the beginning of sequence
const = 10; % interval size of constant pulses after the individual ones
intervals = 3; % number of intervals with constant flips
maxiter = 50; % maximal iterations for algorithm
B1_file = 'B1maps_test4.mat'; % name of B1maps file
pow_constr = 1; % 0 = no constr; 1 = total power; 2 = channel total power
chpow_factor = 1.5; % scale total CHANNEL power factor
pow_factor = 2.0; % scale total power factor
initial = 3; %focus only on states from t=initial; 3 is first echo
% =========================================================================