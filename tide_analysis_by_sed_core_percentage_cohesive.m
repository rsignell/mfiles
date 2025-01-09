% bed_mobility_by_sed_core_percentage

ncclear, clc
load perc_grainsize_crit_smaller_percentiles box_obs_in mc gsizel p

split = find(gsizel < 0.004e-3,1,'last');

gList = 10:10:90;

tide_frac_over = NaN([size(mc,1) size(gList)]);
tide_frac_all = NaN([size(mc,1) size(gList)]);
obs_crit_perc = NaN([size(mc,1) size(gList)]);
obs_percentiles = NaN([size(mc,1) size(gList)]);

mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';
nc = cfdataset(mp.urlH);
tauwc = nc.variable('tauwc');
tauc_tide = nc.variable('tauc_tide');

for gg = 1:size(mc,1)
    ii = box_obs_in(gg,1);
    jj = box_obs_in(gg,2);
    
    if isnan(ii), continue, end
    disp(['On gg = ' num2str(gg) ' of ' num2str(size(mc,1)) '.'])
    
    tSeries = tauwc.data(:,ii,jj);
    tSeries2 = tauc_tide.data(:,ii,jj);
    tSeries(tSeries == 0) = NaN;
    tFrac = tSeries2./tSeries;
    
    gall = mc(gg,:); gall(isnan(gall)) = 0;
    if all(gall==0), continue, end
    gall = cumsum(gall);
    
    for g = 1:length(gList)
        if gall(split) >= 7.5
            obs_crit_perc(gg,g) = 0.1;
        else
            tt = find(gall >= gList(g),1,'first');
            obs_percentiles(gg,g) = gsizel(tt);
            [~,obs_crit_perc(gg,g)] = pmsoulsby(obs_percentiles(gg,g),1);
        end
        tide_frac_over(gg,g) = nanmean(tFrac(tSeries > obs_crit_perc(gg,g)));
    end
            tide_frac_all(gg) = nanmean(tFrac);
    clear gall Tseries
end
clear nc

clear ii jj gg g tt mc mp gsizel box_obs_in ans clear tauwc tauc_tide
obs_crit_perc = squeeze(obs_crit_perc);
% obs_model_percentiles = squeeze(obs_model_percentiles);
tide_frac_over = squeeze(tide_frac_over);
tide_frac_all = squeeze(tide_frac_all);
obs_percentiles = squeeze(obs_percentiles);
save mobility_frequency_cohesive