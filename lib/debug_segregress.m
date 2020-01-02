function debug_segregress(X, Y, R, fig_id, style)

if nargin < 4
    fig_id = 167;
end

if nargin < 5
    style = 'log';
end

[N, D] = size(X);
d = round(D / 2);
Kz = size(R.Z1, 2);
K = length(R.Alpha);

D_vec = zeros(N, 1);
for i = 1:N
    xd = X(i, 1:d);
    xu = X(i, d+1:end);
    D_vec(i) = norm(xd - xu);
end

maxDistance = max(D_vec);

if Kz <= 3
    Gtest = max(log10(Y));
    Ytest = R.Alpha * Gtest + R.Beta;
    [~, I] = sort(Ytest, 'descend');
    R.Alpha = R.Alpha(I);
    R.Beta = R.Beta(I);
    R.Sigma2(1:K) = R.Sigma2(I);
    R.Z1(:, 1:K) = R.Z1(:, I);
    R.Pi(1:K) = R.Pi(I);

    cmap = [0 0 1
            0 1 0
            1 0 0];
    % This is for demonstration in the WI-UAV paper
    if Kz == 2 && K == 2     % Baseline for K = 2
        cmap = [0 0 1
                1 0 0];
    elseif Kz == 3 && K == 3 % Baseline for K = 3
        cmap = [0 0 1
                0 0.5 0
                1 0 0];
    elseif Kz == 3 && K == 2 % proposed one
        cmap = [0 0 1
                1 0 0
                0 0.5 0];
    end
else
    cmap = jet(Kz);
end

figure(fig_id),
for i = 1:N
    [~, kz] = max(R.Z1(i, :));
    if strcmp(style, 'log')
        plot(log10(D_vec(i)), Y(i), 'o', 'color', cmap(kz, :));
    else
        plot((D_vec(i)), Y(i), 'o', 'color', cmap(kz, :));
    end
    hold on
end

dmax = max(D_vec);
dmin = min(D_vec);
dvec = dmin : (dmax - dmin) / 100: dmax;
for k = 1:K
    if strcmp(style, 'log')
        plot(log10(dvec), R.Alpha(k) * log10(dvec) + R.Beta(k), '--', 'color', cmap(k, :));
    else
        plot((dvec), R.Alpha(k) * log10(dvec) + R.Beta(k), '--', 'color', cmap(k, :));
    end
    tx = maxDistance * 0.6;
    ty = R.Alpha(k) * log10(maxDistance) + R.Beta(k);
    text(tx, ty, sprintf('alpha = %5.1f, beta = %5.1f, sigma = %5.1f', ...
                            R.Alpha(k), R.Beta(k), sqrt(R.Sigma2(k))));
end
hold off
if strcmp(style, 'log')
    xlabel('log distance log(d) (m)');
else
    xlabel('distance d (m)');
end
ylabel('channel gain (dB)');
% title('Classification results of the air-to-ground measurement data');
ylim([-130 -50]);
box on
BoxLineWidth = 2;
set(gca, 'LineWidth', BoxLineWidth);
% set(gca, 'XTick', 0:200:1000)
set(gca, 'YTick', -130:20:-50)

    