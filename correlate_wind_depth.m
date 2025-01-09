%% correlate_wind_depth

ncclear, clc, close all
load(fullfile('C:\Users\sdalyander\Documents\StressAnalysis\WavesOnly', ...
    'analysis3_corr_winds'))
windAll = squeeze(xcor_wind_wave(1,:,:));

nc = ncdataset(fullfile('C:\Users\sdalyander\Documents\COAWST\EC30Day\grids', ...
    'USeast_grd17_psd3.nc'));
h = nc.data('h');
clear nc

h(h<=5) = NaN;  %Land

h = h(:);
windAll = windAll(:);

%%
binSize = 10;
depthBins = [0:binSize:200];

numObs = NaN([length(depthBins)-1 1]);
mmean = NaN([length(depthBins)-1 1]);
mstd = NaN([length(depthBins)-1 1]);
x = NaN([length(depthBins)-1 1]);
for ii = 1:length(depthBins)-1
    isIn = find((h >= depthBins(ii)) & (h < depthBins(ii+1)));
    numObs(ii) = length(~isnan(windAll(isIn)));
    mmean(ii) = nanmean(windAll(isIn));
    mstd(ii) = nanstd(windAll(isIn));
    x(ii) = mean(depthBins(ii:ii+1));
end

%%
figure, hold on
for ii = 1:length(x)
    l = line([x(ii) x(ii)], [mmean(ii)-mstd(ii) mmean(ii)+mstd(ii)]);
    l2 = line([depthBins(ii) depthBins(ii+1)], [mmean(ii) mmean(ii)]);
    set(l, 'color', 'r')
    set(l2,'color','b')
    plot(x(ii),mmean(ii),'b.')
end

l = line([100 100], ylim);
set(l,'color','k','linestyle', ':')
box on
xlim([0 200])
xlabel('Depth (m)')
ylabel('r (unitless)')
title('Zero Lag Correlation of Wave Stress to Local Wind Stress')

%%
print('-dpng',fullfile('C:\Users\sdalyander\Documents\StressAnalysis\WavesOnly\figures', ...
    'depth_relationship_corr_local_wind'))
print('-depsc2', fullfile('C:\Users\sdalyander\Documents\StressAnalysis\WavesOnly\figures', ...
    'depth_relationship_corr_local_wind'))