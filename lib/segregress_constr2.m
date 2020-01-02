function R = segregress_constr2(X, Y, K, R0)
% Segmented Regression
% 
% Revision History:
%   Version 5.1 (from segregress_constr.m, K=2 segment version): we
%   modify the intialization method, where we first do one-line linear
%   regression, and then classify the data into two groups according to
%   whether the data is above or below the line. In the EM iteration, there
%   are two phases. In the first phase, the line regression does not have
%   any constraint. In the second phase after the convergence of the first
%   phase, we impose the line constraints. 
%
%   Version 5 (from segregress4.m, constrained version): we add two set of
%   constraints: alpha_k+1 < alpha_k, alpha_k+1 + beta_k+1 < alpha_k +
%   beta_k. The parameters are solved by quadprog.m
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

if K~= 2
    error('This function only supports K=2');
end

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
    
    % Initial clustering via Grand-Regression
    % Step 1: Grand regression
    A = zeros(2);
    b = zeros(2, 1);
    A(1, 1) = sum(G2);
    A(1, 2) = sum(G);
    A(2, 1) = A(1, 2);
    A(2, 2) = N;
    b(1) = sum(G .* Y);
    b(2) = sum(Y);
    
    Lb = [-45, -Inf];
    Ub = [-15, Inf];
    
    options = optimoptions('quadprog',...
                        'Algorithm','interior-point-convex','Display','off');
    Alpha_beta = quadprog(A, -b, [], [], [],[], Lb, Ub, [], options);
    alpha0 = Alpha_beta(1); beta0 = Alpha_beta(2);
    
    % Step 2: Clustering
    Z1 = zeros(N, K);
    Z1(Y >= alpha0 * G + beta0, 1) = 1;
    Z1(Y < alpha0 * G + beta0, 2) = 1;

    % Step 4: Initialize the parameters for the regression model
    Pi = zeros(1, K);   % Marginal probability of the cluster labels
    Alpha = zeros(1, K);
    Beta = zeros(1, K);
    Sigma2 = zeros(1, K);
    
    % ---- New Version (constrained regression) ----
    H = zeros(2 * K); f = zeros(2 * K, 1);
    for k = 1:K
        A = zeros(2);
        b = zeros(2, 1);
        A(1, 1) = sum(Z1(:, k) .* G2);
        A(1, 2) = sum(Z1(:, k) .* G);
        A(2, 1) = A(1, 2);
        A(2, 2) = sum(Z1(:, k));
        b(1) = sum(Z1(:, k) .* G .* Y);
        b(2) = sum(Z1(:, k) .* Y);
        
        H((k - 1) * 2 + [1 2], (k - 1) * 2 + [1 2]) = A;
        f((k - 1) * 2 + [1 2]) = - b;

    end
    
    Lb = [-45, -Inf, -45, -Inf];
    Ub = [-15, Inf, -15, Inf];
    
    options = optimoptions('quadprog',...
                        'Algorithm','interior-point-convex','Display','off');
    % Alpha_beta = quadprog(H, f, Aneq, bneq, [],[],[],[],[],options);
    Alpha_beta = quadprog(H, f, [], [], [],[],Lb,Ub,[],options);
    % ---- New Version (constrained regression) ----
    
    
    for k = 1:K
        Pi(k) = length(find(Z1(:, k) > 0)) / N;

%         % --- Old version ---
%         A = zeros(2);
%         b = zeros(2, 1);
%         A(1, 1) = sum(Z1(:, k) .* G2);
%         A(1, 2) = sum(Z1(:, k) .* G);
%         A(2, 1) = A(1, 2);
%         A(2, 2) = sum(Z1(:, k));
%         b(1) = sum(Z1(:, k) .* G .* Y);
%         b(2) = sum(Z1(:, k) .* Y);
%         alpha_beta = A \ b;
%         % --- old version ---

        % ---- New Version (constrained regression) ----
        Alpha(k) = Alpha_beta((k - 1) * 2 + 1);
        Beta(k) = Alpha_beta((k - 1) * 2 + 2);
        % ---- New Version (constrained regression) ----
        Sigma2(k) = sum( Z1(:, k) .* (Y - Alpha(k) * G - Beta(k)).^2 ) ...
                    / sum( Z1(:, k));
    end
    
%     figure(538), 
%     plot(G(Z1(:, 1) > Z1(:, 2)), Y(Z1(:, 1) > Z1(:, 2)), 'ro'); hold on
%     plot(G(Z1(:, 1) < Z1(:, 2)), Y(Z1(:, 1) < Z1(:, 2)), 'bo'); 
%     Gmin = min(G(:)); Gmax = max(G(:)); Gvec = Gmin:0.05:Gmax;
%     plot(Gvec, alpha0 * Gvec + beta0, 'g--');
%     plot(Gvec, Alpha(1) * Gvec + Beta(1), 'r-.');
%     plot(Gvec, Alpha(2) * Gvec + Beta(2), 'b-.');hold off
    
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
while (norm(Z0 - Z1, 'fro') / norm(Z1, 'fro') > eps || i <= 5) && i < MAXLOOP
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
    
    % ---- New Version (constrained regression) ----
    H = zeros(2 * K); f = zeros(2 * K, 1);
    Aneq = zeros(2 * (K - 1), 2 * K); bneq = zeros(2 * (K - 1), 1);
    for k = 1:K
        A = zeros(2);
        b = zeros(2, 1);
        A(1, 1) = sum(Z1(:, k) .* G2);
        A(1, 2) = sum(Z1(:, k) .* G);
        A(2, 1) = A(1, 2);
        A(2, 2) = sum(Z1(:, k));
        b(1) = sum(Z1(:, k) .* G .* Y);
        b(2) = sum(Z1(:, k) .* Y);
        
        H((k - 1) * 2 + [1 2], (k - 1) * 2 + [1 2]) = A;
        f((k - 1) * 2 + [1 2]) = - b;

        if k < K
            Aneq((k - 1) * 2 + 1, (k - 1) * 2 + (1:4)) = [-1 0 1 0]; % -alpha_k + alpha_k+1 < 0
            Aneq((k - 1) * 2 + 2, (k - 1) * 2 + (1:4)) = [-1.5 -1 1.5 1];% -alpha_k - beta_k + alpha_k+1 + beta_k+1 < 0
        end
    end
    
    Lb = [-45, -Inf, -45, -Inf];
    Ub = [-15, Inf, -15, Inf];
    
    options = optimoptions('quadprog',...
                        'Algorithm','interior-point-convex','Display','off');
    if i >= 5
        Alpha_beta = quadprog(H, f, Aneq, bneq, [],[],Lb,Ub,[],options);
    else
        Alpha_beta = quadprog(H, f, [], [], [],[],Lb,Ub,[],options);
    end
    % ---- New Version (constrained regression) ----
    
    for k = 1:K
        Pi(k) = sum(Z1(:, k)) / N;

%         % --- Old version ---
%         A = zeros(2);
%         b = zeros(2, 1);
%         A(1, 1) = sum(Z1(:, k) .* G2);
%         A(1, 2) = sum(Z1(:, k) .* G);
%         A(2, 1) = A(1, 2);
%         A(2, 2) = sum(Z1(:, k));
%         b(1) = sum(Z1(:, k) .* G .* Y);
%         b(2) = sum(Z1(:, k) .* Y);
%         alpha_beta = A \ b;
%         % --- old version ---

        % ---- New Version (constrained regression) ----
        Alpha(k) = Alpha_beta((k - 1) * 2 + 1);
        Beta(k) = Alpha_beta((k - 1) * 2 + 2);
        % ---- New Version (constrained regression) ----
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


