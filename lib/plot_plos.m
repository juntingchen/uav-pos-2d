DATA = load('../uavDataset/losStatistics.mat');
losStat.Plos = DATA.Plos;
losStat.ElvAngles = DATA.ElvAngles;
clear DATA

figure(1),
plot(losStat.ElvAngles/pi*180, losStat.Plos)
xlabel('Elevation angle [deg]');
ylabel('Empirical LOS probability');
set(gca, 'YTick', 0.2:0.2:1);
set(gca, 'XTick', 0:15:90);
tune_figure;
legend off

