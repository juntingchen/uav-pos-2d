% Massive simulation

clear
addpath(genpath('lib')),

Nue = 10000; % <- reduce this number to shorter simulation time (coarser results)

DATA = load('citymap/urbanMapSingleUserK2.mat');
U = DATA.U; PosBS = DATA.PosBS; 

DATA = load('citymap/losStatistics.mat');
losStat.Plos = DATA.Plos;
losStat.ElvAngles = DATA.ElvAngles;
clear DATA

load('citymap/topologyK2.mat');

U.K = 2;
if U.K == 2
    U.Alpha = [-21.4, -30.3];
    U.Beta =[-36.92, -38.42];
elseif U.K == 3
    U.Alpha = [-22, -28, -36];
    U.Beta =[-28, -24, -22];
else
    error('K should be 2 or 3.');
end
U.A0 = -20.8; U.B0 = -38.5;
U.A1 = U.Alpha(1); U.B1 = U.Beta(1); 
U.A2 = U.Alpha(2); U.B2 = U.Beta(2);


Noise_dBm = -80;
Power_BS_dBm = 33;
Power_UAV_dBm = 33;

U.Noise = 10^(Noise_dBm/10) / 1000; % Watt in linear scale
U.Pb = 10^(Power_BS_dBm/10) / 1000; 
U.Pd = 10^(Power_UAV_dBm/10) / 1000; 

U.Hbs = 45;     % meter, BS height
U.Hmin = 45;    % meter, minimum UAV operation height
U.Hdrone = 50;  % meter, UAV search height
stepSizeMeter = 5;  % UAV search step size

% fun = @(x,y) max(-log2(1 + U.Pd * real(x)), -log2(1 + U.Pb * real(y)));
% fun0 = @(x) -log2(1 + U.Pb * x);

fun = @(x,y) log(1 + 1/(U.Pd * x) + 1/(U.Pb * y));
fun0 = @(x) log(1 + 1/(U.Pb * x));


%%
N_scheme = 6;

tic

Nue = min(size(Topology, 1), Nue);
Rates0 = zeros(Nue, N_scheme);
strongUserIds = zeros(Nue, 1);
failIds = zeros(Nue, 1);

ShuffledIds = randperm(size(Topology, 1));
Topology = Topology(ShuffledIds, :);
parfor i = 1:Nue
    
    PosUE = Topology{i}.PosUE; 
    Blds = Topology{i}.Blds;
    BldTypes = Topology{i}.BldTypes;
    BldLines = Topology{i}.BldLines; 
    BldHeight = Topology{i}.BldHeight; 
    Nbld = size(Blds, 1);
    
    los = IsLosK(PosUE, [PosBS, U.Hbs], BldLines, BldHeight, U.Hdrone, BldTypes);
    if los == 1
        strongUserIds(i) = 1;
        % continue    % We are only interested in the case where the direct BS-user link is blocked
    end
    
    urbanMap = struct();
    urbanMap.BldLines = BldLines;
    urbanMap.BldHeight = BldHeight;
    urbanMap.BldTypes = BldTypes;
    try
        % [Fmin3, Xhat3] = finduavpos3d(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
        Fmin3 = 0;
        [Fmin2, Xhat2] = finduavpos(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
        [Fmin1, Xhat1] = finduavpos1d(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
        % [Fmin_exhst, Xhat_exhst] = finduavpos2d_exhst(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
        Fmin_exhst = Fmin3;
        [FminStat, XhatStat] = finduavposStat(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap, losStat);
    catch
        Fmin1 = 0;
        Fmin2 = 0;
        Fmin3 = 0;
        Fmin_exhst = 0;
        FminStat = 0;
        failIds(i) = 1;
    end
    
    % Direct BS-user link
    k = round((1 - los) * (U.K - 1) + 1);   % propagation segment index
    d = norm([PosBS, U.Hbs] - [PosUE, 0], 2);
    snr = 10 ^ ((U.Alpha(k) * log10(d) + U.Beta(k)) / 10) / U.Noise;
    F0 = fun0(snr);
    
    Rates0(i, :) =  [F0, FminStat, Fmin1, Fmin2, Fmin3, Fmin_exhst];
end

toc

%% Plot results
my_line_styles = {'-', '--', '-.', ':'}.';
Alg_scheme_name = {
    'Direct BS-User link'
    'Probabilistic'
    '1D Optimizationxxx'
    '2D Optimization'
    '3D Optimization'
    'Exhaustive'
};

schemes_to_show = [1 2 3 4 6];

validUserId = failIds < 1;
Rates = Rates0(validUserId, :);

Nue = size(Rates, 1);
n_schemes_to_show = length(schemes_to_show);

X_data = zeros(100, n_schemes_to_show);
F_data = zeros(100, n_schemes_to_show);

Npoints = 1000;
eX_data = zeros(Npoints, n_schemes_to_show);
eF_data = zeros(Npoints, n_schemes_to_show);

for i = 1:n_schemes_to_show
    n = schemes_to_show(i);
    [eF, eX] = ecdf((Rates(:, n) ./ Rates(:, 1)));
    Nef = length(eF);
    I = 1 : (Nef - 1) / Npoints : Nef - (Nef - 1) / (2 * Npoints);
    I = round(I);
    
    eX_data(:, i) = eX(I(:));
    eF_data(:, i) = eF(I(:));
end

h = plot(eX_data(:, 1:end), eF_data(:, 1:end),'linewidth',2);
xlabel('BER(UAV)/BER(BS-USER)');
ylabel('CDF');
set(gca, 'YTick', 0:0.2:1);
set(gca, 'xscale', 'log');
legend(Alg_scheme_name{schemes_to_show(1:end)}, 'location', 'northwest');
xlim([1e-3, 1e-0]);
tune_figure,
% set(h(1), 'linewidth', 2);
% set(h(1), 'Marker', '*', 'Markersize', 6);
% set(h(1), 'LineStyle', ':');
set(h(2), 'LineStyle', '-.');
set(h(3), 'LineStyle', ':');
set(h(4), 'LineStyle', '-');
set(h(4), 'LineWidth', 3);
set(h(5), 'LineStyle', '--');
set(h(5), 'LineWidth', 3);
children = get(gca, 'children');
delete(children(5));
