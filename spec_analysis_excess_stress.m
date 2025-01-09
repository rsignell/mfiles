% bed_mobility_by_sed_core_percentage

ncclear, clc
load perc_grainsize_crit_smaller_percentiles box_obs_in mc gsizel p

split = find(gsizel < 0.004e-3,1,'last');

gList = 10:10:90;

spec_analysis_high = NaN([size(mc,1) size(gList)]);
spec_analysis_low = NaN([size(mc,1) size(gList)]);
spec_analysis_mean = NaN([size(mc,1) size(gList)]);
obs_crit_perc = NaN([size(mc,1) size(gList)]);
obs_percentiles = NaN([size(mc,1) size(gList)]);

mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';
nc = cfdataset(mp.urlH);
tauwc = nc.variable('tauwc');

%time cut off for spectral analysis
Tcut = 33;  %hours
% Hs = spectrum.periodogram;
Hs = spectrum.welch;
for gg = 1:size(mc,1)
    ii = box_obs_in(gg,1);
    jj = box_obs_in(gg,2);
    
    if isnan(ii), continue, end
    disp(['On gg = ' num2str(gg) ' of ' num2str(size(mc,1)) '.'])
    
    total = tauwc.data(:,ii,jj);
    total = R(gg)*total;    %Skin friction only
    
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
        newStress = total - obs_crit_perc(gg,g);
        newStress(newStress < 0) = 0;
        spec_analysis_mean(gg,g) = nanmean(newStress)^2;
        newStress = newStress - nanmean(newStress);
        try
            hmss = msspectrum(Hs,newStress,'Fs',1);
        catch
            continue
        end
        f = hmss.Frequencies;
        pow = hmss.Data;
        low = (f<=1/Tcut);
        high = (f>1/Tcut);
        highfrac = trapz(f(high),pow(high));
        lowfrac = trapz(f(low),pow(low));
        spec_analysis_high(gg,g) = highfrac;
        spec_analysis_low(gg,g) = lowfrac;
        
    end
    clear gall total tide wave curr res
end
clear nc

clear ii jj gg g tt mc mp gsizel box_obs_in ans clear tauwc tauc_tide
obs_crit_perc = squeeze(obs_crit_perc);
spec_analysis_high2mean = spec_analysis_high./spec_analysis_mean;
spec_analysis_low2mean = spec_analysis_low./spec_analysis_mean;
spec_analysis_low2high = spec_analysis_low./spec_analysis_high;
clear tauc tauc_res tauw split highfrac lowfrac low high f hmss newStress pow

save mobility_frequency