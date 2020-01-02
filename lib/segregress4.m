function R = segregress4(X, Y, K, R0)
% Segmented Regression
% 
% Revision History:
%   Version 2 (changed name to segreress4.m): Develop new method for 
%     Globecom Wi-UAV workshop paper revision. New ideas:
%       - Clustering algorithm, added constraints to the clustering, only
%       cluster to consistent segments
%       - K + 1 semgents, where the last segment is to incorporate the
%       background noise (Two types of noise: segment noise and background
%       noise)
%   Version 1.1: Added initial value R0; ordered the segments
%   Version 1: Segmented regression function for ICC'17 paper
%
% Segmented regression for the channel model based on maximul likelihood
% (ML) estimation and expectation-maximizaiton (EM) algorithm
%
% INPUT
%   X           N * D matrix, for N data samples; each row contains a
%               sample with d dimension. 
%   Y           N * 1 vector, for N data samples. 
%   K           The number of segments/clusters/regions
%
% OUTPUT
%   R.Alpha
%   R.Beta
%   R.Sigma2
%   R.Pi
%   R.Xc

[N, D] = size(X);
G = zeros(N, 1);
d = round(D / 2);
for l = 1:N
    G(l) = log10(norm(X(l, 1:d) - X(l, d+1:end)));
end
G2 = G .^ 2;
    
var_background_noise = 100;
k_bg = ceil(K / 2);

% Initialization ----------------------------------------------------------
if nargin >= 4
    % The case when the initial value is provided 
    Alpha = [R0.Alpha(1:K) R0.alpha_background];
    Beta = [R0.Beta(1:K) R0.beta_background];
    Sigma2 = [R0.Sigma2(1:K) R0.sigma2_background];
    Pi = R0.Pi;
    
    % Initial classification based on new dataset (X, Y)
    Z1 = zeros(N, K + 1);
    Pxy = zeros(N, K + 1);      % The joint density p_k(x, y) for the k-th cluster eval @ l-th data sample (x, y)
    for l = 1:N
        for k = 1:K + 1
            Pxy(l, k) = exp( - (Y(l) - Alpha(k) * G(l) - Beta(k))^2 / 2 / Sigma2(k)) ...
                        / (sqrt(2 * pi * Sigma2(k)));
        end
        for k = 1:K + 1
            Z1(l, k) = Pxy(l, k) * Pi(k) / sum(Pxy(l, :) .* Pi);
        end
    end
    
else
    
    % Initial clustering via Kmeans 
    % Step 1: Regularize the data (for same variance of X and Y)
    varY = var(Y);
    varX = var(X(:, 1));
    Y1 = Y * sqrt(varX / varY);
    W = [X Y1];
    Z = zeros(N, K);    % cluster label

    % Step 2: Initialization
    % Randomly assign K data samples into K clusters
    I = randperm(N);
    Wc = zeros(K, D + 1);
    for k = 1:K
        Z(I(k), k) = 1;            % I(k)-th data sample to represent the k-th cluster
        Wc(k, :) = W(k, :);     % Cluster center
    end

    % Step 3: K-means Iteration
    MAXLOOP = 10;
    i = 0;
    Z1 = Z;
    Z0 = zeros(N, K);
    while norm(Z0(:) - Z1(:)) > 0 && i < MAXLOOP
        i = i + 1;

        Z0 = Z1;
        Z1 = zeros(N, K);
        % Nearest neighbor assignment
        for l = 1:N
            % For each data sample, assign it to the cluster with the shortest
            % distance to the cluster center
            distances = zeros(1, K);
            for k = 1:K
                distances(k) = norm(W(l, :) - Wc(k, :));
            end
            [~, I] = min(distances);
            Z1(l, I) = 1;
        end

        % Update the cluster center
        for k = 1:K
            I = find(Z1(:, k) > 0);
            if ~isempty(I)
                Wc(k, :) = mean(W(I, :), 1);
            end
        end
    end

    % Step 4: Initialize the parameters for the regression model
    Pi = zeros(1, K);   % Marginal probability of the cluster labels
    Alpha = zeros(1, K);
    Beta = zeros(1, K);
    Sigma2 = zeros(1, K);
    for k = 1:K
        Pi(k) = length(find(Z1(:, k) > 0)) / N;

        A = zeros(2);
        b = zeros(2, 1);
        A(1, 1) = sum(Z1(:, k) .* G2);
        A(1, 2) = sum(Z1(:, k) .* G);
        A(2, 1) = A(1, 2);
        A(2, 2) = sum(Z1(:, k));
        b(1) = sum(Z1(:, k) .* G .* Y);
        b(2) = sum(Z1(:, k) .* Y);

        alpha_beta = A \ b;

        Alpha(k) = alpha_beta(1);
        Beta(k) = alpha_beta(2);
        Sigma2(k) = sum( Z1(:, k) .* (Y - Alpha(k) * G - Beta(k)).^2 ) ...
                    / sum( Z1(:, k));
    end
    
end

% EM algorithm ------------------------------------------------------------
% We create a background segmente labeled as K + 1. The background is
% modeled using the ceili(K / 2) segment, with a large enough and fixed 
% variance (here we set it as 20 dB (std). 
Z0 = zeros(N, K + 1);
if nargin >= 4
else 
    Alpha = [Alpha, Alpha(k_bg)];
    Beta = [Beta, Beta(k_bg)];
    Sigma2 = [Sigma2, var_background_noise];
    Z1 = [Z1(:, 1:K), zeros(size(Z1, 1), 1)];
    Pi = [Pi(1:K) * 0.95, 0.05];
end
eps = 1e-3;
i = 0;

MAXLOOP = 200;
Q = zeros(1, MAXLOOP);
while norm(Z0 - Z1, 'fro') / norm(Z1, 'fro') > eps && i < MAXLOOP
    i = i + 1;
    
    % E-step: Update the soft label (i.e., the probably Z(l, k) that the l-th
    % sample is assigned to the k-th cluster
    Z0 = Z1;
    Z1 = zeros(N, K + 1);
    Pxy = zeros(N, K + 1);      % The joint density p_k(x, y) for the k-th cluster eval @ l-th data sample (x, y)
    for l = 1:N
        for k = 1:K + 1
            Pxy(l, k) = exp( - (Y(l) - Alpha(k) * G(l) - Beta(k))^2 / 2 / Sigma2(k)) ...
                        / (sqrt(2 * pi * Sigma2(k)));
        end
        for k = 1:K + 1
            Z1(l, k) = Pxy(l, k) * Pi(k) / sum(Pxy(l, :) .* Pi);
        end
    end

    % M-step: Update label marginal probabiliy and the regresssion parameters
    for k = 1:K
        Pi(k) = sum(Z1(:, k)) / N;

        A = zeros(2);
        b = zeros(2, 1);
        A(1, 1) = sum(Z1(:, k) .* G2);
        A(1, 2) = sum(Z1(:, k) .* G);
        A(2, 1) = A(1, 2);
        A(2, 2) = sum(Z1(:, k));
        b(1) = sum(Z1(:, k) .* G .* Y);
        b(2) = sum(Z1(:, k) .* Y);

        alpha_beta = A \ b;

        Alpha(k) = alpha_beta(1);
        Beta(k) = alpha_beta(2);
        Sigma2(k) = sum( Z1(:, k) .* (Y - Alpha(k) * G - Beta(k)).^2 ) ...
                    / sum( Z1(:, k));
    end

    % Background segment
    Pi(K + 1) = sum(Z1(:, K + 1)) / N;
    Alpha(K + 1) = Alpha(k_bg);
    Beta(K + 1) = Beta(k_bg);
    Sigma2(K + 1) = var_background_noise;
    
    Q(i) = sum(sum(Z1 .* (repmat(log(Pi + 1e-11), N, 1) + log(Pxy + 1e-11))));
end 

% Arrange the label indices
Gtest = max(G);
alpha_background = Alpha(K + 1);
beta_background = Beta(K + 1);
sigma2_background = Sigma2(K + 1);

Ytest = Alpha(1:K) * Gtest + Beta(1:K);
[~, I] = sort(Ytest, 'descend');
Alpha = Alpha(I);
Beta = Beta(I) ;
Sigma2 = Sigma2(I);
Z1 = Z1(:, [I, K + 1]);
Pi = Pi([I, K + 1]);

R.Alpha = Alpha;
R.Beta = Beta;
R.Sigma2 = Sigma2;
R.alpha_background = alpha_background;
R.beta_background = beta_background;
R.sigma2_background = sigma2_background;
R.Pi = Pi;
R.Q = Q(1:i);
R.Z1 = Z1;
R.X = X;


