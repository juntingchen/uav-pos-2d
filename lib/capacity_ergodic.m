function R = capacity_ergodic(Ks_dB, SNRs_dB)
% Generate the Ergodic capacity as a function of average SNR for fading
% channels.
%
% INPUT: 
%   Ks      K-factors. K=-Inf means Rayleigh fading. # of cases = length(Ks)
%   SNRs    A sequence of SNR values
%
% OUTPUT:
%   Rawgn   AWGN channel capacity evaluated over the sequence of SNR values
%   R       length(K) * length(SNRs) matrix, returning the Ergodic capacity
%           evaluated using Monte-Carlo simulations

if nargin < 1
    Ks_dB = [-Inf, 9, 20];
    SNRs_dB = [-10:2:20];
    str_legend = {'AWGN', 'Rayleigh', '9 dB', '20 dB'};
end


Nmc = 10000;
M = length(Ks_dB);
N = length(SNRs_dB);

Rawgn = zeros(1, N);
Rerg = zeros(M, N);

for j = 1:N
    snr = 10^(SNRs_dB(j)/10);
    Rawgn(j) = log2(1 + snr);
end
    
for i = 1:M
    K = 10^(Ks_dB(i)/10);
    if K < 1e-10
        K = 1e-10;
    end
    
    sigma = sqrt(1/(2 * (1+K)));
    nu = sqrt(K / (K + 1));
    
    amplitude_gains = random('Rician', nu, sigma, [1, Nmc]);
    
    for j = 1:N
        snr = 10^(SNRs_dB(j)/10);
        Rerg(i, j) = mean(log2(1 + snr * amplitude_gains.^2));
    end
end

R = [Rawgn
     Rerg];
 
if nargin < 1
    figure(1),
    plot(SNRs_dB, R, 'linewidth', 2);
    legend(str_legend);
    grid on

    hold on
    rate = @(x)ppval(spline(SNRs_dB, Rerg(1, :)), x);
    x_snr = -20:30;
    plot(x_snr, rate(x_snr));
    hold off
end
