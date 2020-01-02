% Massive simulation

clear
addpath(genpath('lib')),

Nue = 10000;    % <- reduce this number to shorter simulation time (coarser results)

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

fun = @(x,y) max(-log2(1 + U.Pd * real(x)), -log2(1 + U.Pb * real(y)));
fun0 = @(x) -log2(1 + U.Pb * x);

% Ergodic capacity
SNRs_dB = -10:2:20; Ks_dB = [9, -Inf];
Rerg = capacity_ergodic(Ks_dB, SNRs_dB);

fun1 = @(x,y) max(- max(0, ppval(spline(SNRs_dB, Rerg(1, :)), 10 * log10(U.Pd * real(x)))), ... 
                  - log2(1 + U.Pb * real(y))); % UAV-UE_LOS(K-factor = 9dB, 
              
fun2 = @(x,y) max(- max(0, ppval(spline(SNRs_dB, Rerg(2, :)), 10 * log10(U.Pd * real(x)))), ... 
                  - log2(1 + U.Pb * real(y))); % UAV-UE_NLOS, Rayleigh fading 

%%
N_scheme = 6;

tic

Nue = min(size(Topology, 1), Nue);
Rates0 = zeros(Nue, N_scheme);
strongUserIds = zeros(Nue, 1);
failIds = zeros(Nue, 1);
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
    
    % Direct BS-user link
    k = round((1 - los) * (U.K - 1) + 1);   % propagation segment index
    d = norm([PosBS, U.Hbs] - [PosUE, 0], 2);
    snr = 10 ^ ((U.Alpha(k) * log10(d) + U.Beta(k)) / 10) / U.Noise;
    F0 = fun0(snr);
    
    try
        % [Fmin3, Xhat3] = finduavpos3d(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
        % Fmin3 = min(Fmin3, F0);
        Fmin3 = 0;
        
        [~, Xhat2] = finduavpos(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
        los = IsLosK(PosUE, [Xhat2, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
        Fmin2 = getcostf2DK_ergodic([Xhat2, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun1, fun2);
        % Fmin2 = min(Fmin2, F0);
        
        [~, Xhat1] = finduavpos1d(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
        los = IsLosK(PosUE, [Xhat1, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
        Fmin1 = getcostf2DK_ergodic([Xhat1, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun1, fun2);

        % Fmin1 = min(Fmin1, F0);
        
        % [Fmin_exhst, Xhat_exhst] = finduavpos2d_exhst(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
        Fmin_exhst = Fmin3;
        
        [~, XhatStat] = finduavposStat(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap, losStat);
        los = IsLosK(PosUE, [XhatStat, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
        FminStat = getcostf2DK_ergodic([XhatStat, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun1, fun2);
        % FminStat = min(FminStat, F0);
    catch
        Fmin1 = 0;
        Fmin2 = 0;
        Fmin3 = 0;
        Fmin_exhst = 0;
        FminStat = 0;
        failIds(i) = 1;
    end
   
    Rates0(i, :) = - [F0, FminStat, Fmin1, Fmin2, Fmin3, Fmin_exhst];
end

toc

%% Plot results
my_line_styles = {'-', '--', '-.', ':'}.';
Alg_scheme_name = {
    'Direct BS-User linkxx'
    'Probabilistic Alg'
    'Simple Search'
    'Proposed'
    'Proposed (3D)'
    'Exhaustive'
};

schemes_to_show = [1 2 3 4 6];
N_scheme_to_show = length(schemes_to_show);

validUserId = failIds < 1;
Rates = Rates0(validUserId, :);

Nue = size(Rates, 1);

maxdata = max(Rates(:));
Npt = 40;
XI = sort([0.1 0.17 0.3 0.5 (0:1/(Npt - 1 - 4):1) * maxdata], 'ascend');

X_data = zeros(Npt, N_scheme_to_show);
F_data = zeros(Npt, N_scheme_to_show);

for i = 1:N_scheme_to_show
    n = schemes_to_show(i);
    
    r_vec = Rates(:, n);

    [F1,X1] = ksdensity(r_vec, XI, 'function', 'cdf');
    
    X_data(:, i) = X1(:);
    F_data(:, i) = F1;
    
end

figure(1),
h = plot(X_data, F_data,'linewidth', 2);
set(gca, 'FontSize', 14);
legend(Alg_scheme_name{schemes_to_show}, 'location', 'southeast');
xlim([0 ceil(max(Rates(:)))]);
set(gca, 'YTick', 0:0.2:1);
xlabel('bps/Hz');
ylabel('CDF');
tune_figure,
set(h(1), 'linewidth', 2);
set(h(1), 'Marker', '*', 'Markersize', 6);
set(h(1), 'LineStyle', ':');
set(h(2), 'LineStyle', '-.');
set(h(3), 'LineStyle', ':');
set(h(4), 'LineStyle', '-');
set(h(4), 'LineWidth', 3);
set(h(5), 'linestyle', '--');
set(h(5), 'LineWidth', 3);

% ----
schemes_to_show = [1 2 3 4];

figure(2),
rateNoUav = Rates(:, 1);
[~, sortedIndex] = sort(rateNoUav, 'ascend');
low20percentileIndex = sortedIndex(1:round(Nue * 0.2));
high20percentileIndex = sortedIndex(round(Nue * 0.8): end);

RateLow = mean(Rates(low20percentileIndex, schemes_to_show), 1);
RateMean = mean(Rates(:, schemes_to_show), 1);
RateHigh = mean(Rates(high20percentileIndex, schemes_to_show), 1);
h = bar([RateLow
         RateMean
         RateHigh]);
set(gca, 'FontSize', 14);
set(h, 'linewidth', 2);
ylim([0, 10]);
legend(Alg_scheme_name{schemes_to_show}, 'location', 'northwest');
set(gca, 'XTickLabel', {'20th percentile', 'Mean', 'Top 20th percentile'});
set(gca, 'YTick', 0:2:10);
ylabel('Average end-to-end throughput [bps/Hz]');
% label the bars
Xdata = [RateLow
         RateMean
         RateHigh];
bartext = [];
for i = 1:size(Xdata, 1)
    for j = 1:size(Xdata, 2)
        bartext(i, j) = text(i + (j - 2.5) * 0.18, Xdata(i, j) + 0.05, ...
            sprintf('%1.2f', Xdata(i, j)), 'fontsize', 12);
    end
end
% Use the handles TH to modify some properties
set(bartext,'Horizontalalignment','center',...
'verticalalignment','bottom') ;
tune_figure,
[im_hatch,colorlist] = applyhatch_pluscolor(gcf,'\-x./+',0,0,[],150,2,2);

