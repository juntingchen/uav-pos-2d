% POSSIBLY OUTDATED DESCRIPTION
% Multiuser case
% Massive simulation output
%   - rate versus radius
%   - rate versus user group --> Rate: center user location * radius
%
% Parameters:
%   - Radius of a user cluster
%   - user density (e.g., at most 20 users within the radius?)
%
% Figures:
%   - one figure for strong user: Mean rate versus radius
%   - one figure for weak user: Mean rate versus radius
%
% Data structure:
%   -- Rates0_cell ---
%   {cell - 1} center user location 1
%       [cluster radius 1: scheme-1, scheme-2, scheme-3, ...]
%       [cluster radius 2: scheme-1, scheme-2, scheme-3, ...]
%       ...
%       [cluster radius N_radius: scheme-1, scheme-2, scheme-3, ...]
%   {Cell - 2} center user location 2
%   ...
%
clear
addpath(genpath('lib')),

Nue = 100;

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

Radius_vec = [0:5:50];
N_radius = length(Radius_vec);
Nuser_per_cluster = 10;

U.Hbs = 45;     % meter, BS height
U.Hmin = 45;    % meter, minimum UAV operation height
U.Hdrone = 50;  % meter, UAV search height
stepSizeMeter = 5;  % UAV search step size

% % End-to-end data rate
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

% Preprocessing of user topology
maxNumUsers = length(Topology);
AllUEPos = zeros(maxNumUsers, 2);
for i = 1:maxNumUsers
    AllUEPos(i, :) = Topology{i}.PosUE; 
end

tic

Nue = min(size(Topology, 1), Nue);
n_chunk = 10;
Nue_per_chunk = floor(Nue / n_chunk);
Nue = Nue_per_chunk * n_chunk;



Rates_cell = cell(Nue, 1);
strongUserIds = cell(n_chunk, 1);
failIds = cell(n_chunk, 1);

% try
%     poolobj = parpool(4);
% catch
%     warning('fail to initialize parpool!');
% end


for outer_i = 1:n_chunk
    fprintf('Running %d/%d ...\n', outer_i, n_chunk);
    
    Rates0_cell = cell(Nue_per_chunk, 1);
    strongUserIds0 = zeros(Nue_per_chunk, 1);
    failIds0 = zeros(Nue_per_chunk, 1);
    
    parfor i = 1:Nue_per_chunk

        % Cluster central user (representative user)
        PosUE = Topology{i + (outer_i - 1) * Nue_per_chunk}.PosUE; 
        Blds = Topology{i + (outer_i - 1) * Nue_per_chunk}.Blds;
        BldTypes = Topology{i + (outer_i - 1) * Nue_per_chunk}.BldTypes;
        BldLines = Topology{i + (outer_i - 1) * Nue_per_chunk}.BldLines; 
        BldHeight = Topology{i + (outer_i - 1) * Nue_per_chunk}.BldHeight; 
        Nbld = size(Blds, 1);

        los = IsLosK(PosUE, [PosBS, U.Hbs], BldLines, BldHeight, U.Hdrone, BldTypes);
        if los == 1
            strongUserIds0(i) = 1;
            % continue    % We are in particular interested in the case where the direct BS-user link is blocked
        end   

        urbanMap = struct();
        urbanMap.BldLines = BldLines;
        urbanMap.BldHeight = BldHeight;
        urbanMap.BldTypes = BldTypes;
        try
            % [Fmin3, Xhat3] = finduavpos3d(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
            Fmin3 = 0;
            % Fmin3 = 0;
            
            [~, Xhat2, ~, Xhat2Array] = finduavpos(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
            los = IsLosK(PosUE, [Xhat2, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
            Fmin2 = getcostf2DK_ergodic([Xhat2, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun1, fun2);
            
            [~, Xhat1] = finduavpos1d(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
            los = IsLosK(PosUE, [Xhat1, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
            Fmin1 = getcostf2DK_ergodic([Xhat1, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun1, fun2);

            % [Fmin_exhst, Xhat_exhst] = finduavpos2d_exhst(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap);
            Fmin_exhst = Fmin3;
            
            [~, XhatStat] = finduavposStat(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap, losStat);
            los = IsLosK(PosUE, [XhatStat, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
            FminStat = getcostf2DK_ergodic([XhatStat, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun1, fun2);
        

        catch
            Fmin1 = 0; Fmin2 = 0; Fmin3 = 0; Fmin_exhst = 0; FminStat = 0;
            XhatStat = PosUE; Xhat2Array = PosUE;
            failIds0(i) = 1;
            warning('something is wrong!');
        end

        % Direct BS-user link
        k = round((1 - los) * (U.K - 1) + 1);   % propagation segment index
        d = norm([PosBS, U.Hbs] - [PosUE, 0], 2);
        snr = 10 ^ ((U.Alpha(k) * log10(d) + U.Beta(k)) / 10) / U.Noise;
        F0 = fun0(snr);

        Rate_center_user = - [F0, FminStat, Fmin1, Fmin2, Fmin3, Fmin_exhst];
        Rates_rad = zeros(N_radius, N_scheme);
        Rates_rad(1, :) = Rate_center_user * Nuser_per_cluster;

        for ir = 2:N_radius
            radius = Radius_vec(ir);
            Rates_rad(ir, :) = Rate_center_user;

            % Find other users within the cluster
            UserIds = zeros(1, Nuser_per_cluster);
            UserIds(1) = i + (outer_i - 1) * Nue_per_chunk;
            cnt1 = 1;

            cnt = Nuser_per_cluster;
            randomId = randperm(maxNumUsers);
            cluUserId = zeros(1, Nuser_per_cluster - 1);
            cntTrial = min(Nue, maxNumUsers / 2);
            while cnt > 0
                cntTrial = cntTrial + 1;
                if norm(PosUE - Topology{randomId(cntTrial)}.PosUE) < radius
                    cnt1 = cnt1 + 1;
                    UserIds(cnt1) = randomId(cntTrial);

                    PosMember = Topology{randomId(cntTrial)}.PosUE;
                    PosMember3 = [PosMember, 0];
                    PosBS3 = [PosBS, U.Hbs];
                    Blds = Topology{randomId(cntTrial)}.Blds;
                    BldTypes = Topology{randomId(cntTrial)}.BldTypes;
                    BldLines = Topology{randomId(cntTrial)}.BldLines; 
                    BldHeight = Topology{randomId(cntTrial)}.BldHeight; 
                    Nbld = size(Blds, 1);

                    % los1 = IsLosK(PosMember3, [Xhat1, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
                    % F1 = getcostf2DK([Xhat1, U.Hdrone], PosMember3, PosBS3, los1, U, fun);

                    % los2 = IsLosK(PosMember3, [Xhat2, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
                    % F2 = getcostf2DK([Xhat2, U.Hdrone], PosMember3, PosBS3, los2, U, fun);

                    los0 = IsLosK(PosMember3, [XhatStat, U.Hdrone], BldLines, BldHeight, U.Hdrone, BldTypes);
                    Fstat = getcostf2DK([XhatStat, U.Hdrone], PosMember3, PosBS3, los0, U, fun);

                    % Direct BS-user link
                    los = IsLos(PosMember, [PosBS, U.Hbs], BldLines, BldHeight, U.Hdrone);
                    k = round((1 - los) * (U.K - 1) + 1);   % propagation segment index
                    d = norm([PosBS, U.Hbs] - [PosMember, 0], 2);
                    snr = 10 ^ ((U.Alpha(k) * log10(d) + U.Beta(k)) / 10) / U.Noise;
                    F0 = fun0(snr);
                    Rates_rad(ir, :) = Rates_rad(ir, :) + (- [F0, Fstat, 0, 0, 0, 0]);
                    cnt = cnt - 1;
                end
                if cntTrial >= maxNumUsers
                    cntTrial = 0;
                end
            end

            % BEGIN - proposed search trajectory for multiuser case
            nStep = size(Xhat2Array, 1);

            Fmin = inf;
            Xhat2mu = [0, 0, 0];
            for pos_i = 1:nStep
                Xpos = Xhat2Array(pos_i, :);

                % -- compute the sum cost for all users --
                sum_cost = 0;
                for imem = 1:Nuser_per_cluster
                    uid = UserIds(imem);

                    PosMember = Topology{uid}.PosUE;
                    PosMember3 = [PosMember, 0];
                    PosBS3 = [PosBS, U.Hbs];
                    Blds = Topology{uid}.Blds;
                    BldTypes = Topology{uid}.BldTypes;
                    BldLines = Topology{uid}.BldLines; 
                    BldHeight = Topology{uid}.BldHeight; 
                    Nbld = size(Blds, 1);

                    los = IsLosK(PosMember3, Xpos, BldLines, BldHeight, U.Hdrone, BldTypes);
                    f = getcostf2DK([Xpos, U.Hdrone], PosMember3, PosBS3, los, U, fun);
                    sum_cost = sum_cost + f;
                end    

                if sum_cost < Fmin
                    Fmin = sum_cost;
                    Xhat2mu = Xpos;
                end

            end        
            % END - proposed multiuser case
             Rates_rad(ir, 4) = - Fmin;

        end

        Rates0_cell{i} = Rates_rad / Nuser_per_cluster;
    end % parfor i
    
    Rates_cell((outer_i - 1) * Nue_per_chunk + 1:outer_i * Nue_per_chunk) ...
        = Rates0_cell;
    strongUserIds{outer_i} = strongUserIds0;
    failIds{outer_i} = failIds0;
end % outer_i

% delete(poolobj);
toc

failIds = cell2mat(failIds);
strongUserIds = cell2mat(strongUserIds);
Rates0_cell = Rates_cell;

%% Plot results
my_line_styles = {'-', '--', '-.', ':'}.';
my_markers = {'o', 'd', '^', 's', 'v', '+'}.';
Alg_scheme_name = {
    'Direct BS-User link'
    'Probabilistic'
    '1D Optimization'
    '2D Optimization'
    '3D Optimization'
    'Exhaustive'
};

% 
schemes_to_show = [4 2 1];

validStrongUserId = find(strongUserIds >= 1 & failIds < 1);
validWeakUserId = find(strongUserIds < 1 & failIds < 1);

% Strong user group ---------------------------------
RateS = zeros(N_radius, N_scheme);
for i = 1:length(validStrongUserId)
    user_id = validStrongUserId(i);
    RateS = RateS + Rates0_cell{user_id};
end
RateS = RateS / length(validStrongUserId);

figure(1),
h_curve = plot(Radius_vec, RateS(:, schemes_to_show), 'MarkerSize', 9);
set(gca, 'FontSize', 14);
set(h_curve, {'Marker'}, my_markers(1:length(schemes_to_show)));
legend(Alg_scheme_name{schemes_to_show}, 'location', 'northwest');
xlabel('Cluster radius [m]');
ylabel('Average end-to-end throughput [bps/Hz]');
ylim([0 inf]);
% set(gca, 'XTickLabel', {'Weak Users', 'Strong Users'});
% set(gca, 'YTick', 0:1:4);
tune_figure;

% Weak user group ---------------------------------
RateW = zeros(N_radius, N_scheme);
for i = 1:length(validWeakUserId)
    user_id = validWeakUserId(i);
    RateW = RateW + Rates0_cell{user_id};
end
RateW = RateW / length(validWeakUserId);

figure(2),
h_curve = plot(Radius_vec, RateW(:, schemes_to_show));
set(h_curve, {'Marker'}, my_markers(1:length(schemes_to_show)));
set(gca, 'FontSize', 14);
legend(Alg_scheme_name{schemes_to_show}, 'location', 'northeast');
xlabel('Cluster radius [m]');
ylabel('Average end-to-end throughput [bps/Hz]');
maxy = ceil(max(vec(RateW(:, schemes_to_show))) / 1) * 1;
ylim([0 maxy]);
% set(gca, 'XTickLabel', {'Weak Users', 'Strong Users'});
% set(gca, 'YTick', 0:1:4);
tune_figure;

