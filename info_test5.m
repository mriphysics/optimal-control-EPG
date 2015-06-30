% input data for test 5
%% ================ INITIALIZE PARAMETERS =================================
ESP = 12.0 ;%<- echo spacing, ms
T1 = 1500; %<- T1,ms
T2 = 100;  %<- T2,ms
klim = 25; % order of maximum k, use 25
max_a = 1.5; % maximum flip = max_a*pi;
red_b1 = 1; % spatial undersampling factor for B1 maps
a0 = pi/180*[105,174.06, 145.08,140.04,140.04,140.04,...
            140.04,140.04,140.04,140.04]; % target sequence 
% define piecewise constant indices
singles = 10; % number of individual pulses at the beginning of sequence
const = 20; % interval size of constant pulses after the individual ones
intervals = 0; % number of intervals with constant flips
maxiter = 50; % maximal iterations for algorithm
B1_file = 'B1maps_test5.mat'; % name of B1maps file
pow_constr = 1; % 0 = no constr; 1 = total power; 2 = channel total power
chpow_factor = 1.5; % scale total CHANNEL power factor
pow_factor = 2.0; % scale total power factor
initial = 3; %focus only on states from t=initial; 3 is first echo
% =========================================================================