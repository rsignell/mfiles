% bed_mobility_by_sed_core_percentage

load perc_grainsize_crit_smaller_percentiles box_obs_in mc gsizel p

gList = [10:10:90];

obs_percentiles = NaN([size(mc,1) size(gList)]);
obs_crit_perc = NaN([size(mc,1) size(gList)]);
obs_model_percentiles = NaN([size(mc,1) size(gList)]);

mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';
nc = cfdataset(mp.urlH);
var = nc.variable('tauwc');

for gg = 1:size(mc,1)
    ii = box_obs_in(gg,1);
    jj = box_obs_in(gg,2);
    
    if isnan(ii), continue, end
    
    gall = mc(gg,:); gall(isnan(gall)) = 0;
    if all(gall==0), continue, end
    gall = cumsum(gall);
        
    
    tSeries = var.data(:,ii,jj);
    for g = 1:length(gList)
        tt = find(gall >= gList(g),1,'first');
        obs_percentiles(gg,g) = gsizel(tt);
        [~,obs_crit_perc(gg,g)] = pmsoulsby(obs_percentiles(gg,g),1);
           
        obs_model_percentiles(gg,g) = 100*...
            (length(find(tSeries >= obs_crit_perc(gg,g)))./...
            length(find(~isnan(tSeries))));
    end
    clear gall tSeries
end
clear nc

clear ii jj gg g var tt mc mp gsizel box_obs_in ans
obs_crit_perc = squeeze(obs_crit_perc);
obs_model_percentiles = squeeze(obs_model_percentiles);
obs_percentiles = squeeze(obs_percentiles);
save compare_sediment_percentiles_to_stress_tseries