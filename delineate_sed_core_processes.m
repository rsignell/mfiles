% bed_mobility_by_sed_core_percentage

ncclear, clc
load perc_grainsize_crit_smaller_percentiles box_obs_in mc gsizel p

split = find(gsizel < 0.004e-3,1,'last');

gList = 10:10:90;

tide_mean = NaN([size(mc,1) size(gList)]);
wave_mean = NaN([size(mc,1) size(gList)]);
curr_mean = NaN([size(mc,1) size(gList)]);
res_mean = NaN([size(mc,1) size(gList)]);
total_mean = NaN([size(mc,1) size(gList)]);
obs_crit_perc = NaN([size(mc,1) size(gList)]);
obs_percentiles = NaN([size(mc,1) size(gList)]);

mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';
nc = cfdataset(mp.urlH);
tauwc = nc.variable('tauwc');
tauc_tide = nc.variable('tauc_tide');
tauc_res = nc.variable('tauc_res');
tauc = nc.variable('tauc');
tauw = nc.variable('tauw');


for gg = 1:size(mc,1)
    ii = box_obs_in(gg,1);
    jj = box_obs_in(gg,2);
    
    if isnan(ii), continue, end
    disp(['On gg = ' num2str(gg) ' of ' num2str(size(mc,1)) '.'])
    
    total = tauwc.data(:,ii,jj);
    tide = tauc_tide.data(:,ii,jj);
    wave = tauw.data(:,ii,jj);
    curr = tauc.data(:,ii,jj);
    res = tauc_res.data(:,ii,jj);
    
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
        isOver = total > obs_crit_perc(gg,g);
        tide_mean(gg,g) = nanmean(tide(isOver));
        total_mean(gg,g) =  nanmean(total(isOver));
        wave_mean(gg,g) = nanmean(wave(isOver));
        res_mean(gg,g) = nanmean(res(isOver));
        curr_mean(gg,g) = nanmean(curr(isOver));
    end
    clear gall total tide wave curr res
end
clear nc

clear ii jj gg g tt mc mp gsizel box_obs_in ans clear tauwc tauc_tide
obs_crit_perc = squeeze(obs_crit_perc);
tide_mean = squeeze(tide_mean);
res_mean = squeeze(res_mean);
curr_mean= squeeze(curr_mean);
wave_mean = squeeze(wave_mean);
total_mean = squeeze(total_mean);
clear tauc tauc_res tauw split 
save mobility_process