% bed_mobility_by_sed_core_percentage

ncclear, clc
load grainsize_obs box_obs_in mc gsizel p

split = find(gsizel < 0.004e-3,1,'last');

gList = 10:10:90;

mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

obs_percentiles = NaN([size(mc,1) size(gList)]);
obs_crit_perc = NaN([size(mc,1) size(gList)]);
obs_perc_over_pt1 = NaN([size(mc,1) length(mlists)]);
obs_perc_over_pt2 = NaN([size(mc,1) length(mlists)]);
obs_model_percentiles = NaN([size(mc,1) length(mlists) size(gList)]);
R = NaN([size(mc,1) 1]);

mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CompGrid/swan_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\WavesOnly';
nc = cfdataset(mp.urlH);
var = nc.variable('tauw');
time = nc.time('time');
gTime = datevec(time);

for gg = 1:size(mc,1)
    ii = box_obs_in(gg,1);
    jj = box_obs_in(gg,2);
    
    if isnan(ii), continue, end
    disp(['On gg = ' num2str(gg) ' of ' num2str(size(mc,1)) '.'])
    
    tSeries = var.data(:,ii,jj);
    
    gall = mc(gg,:); gall(isnan(gall)) = 0;
    if all(gall==0), continue, end
    gind = gall;
    gall = cumsum(gall);
            
    gind(gall < 90) = NaN;
    d90 = nansum(gind.*gsizel)/nansum(gind);
    if d90 > 0.5e-2
        %Smokies...bigger than the original!
        R(gg) = 1;
    else
    frac = nsol(0.5,@(R)eb52(R,d90,0.5e-2),[eps eps]);
       R(gg) = frac;
    end
    
    for g = 1:length(gList)
        if gall(split) >= 7.5
            obs_crit_perc(gg,g) = 0.1;
        else
            tt = find(gall >= gList(g),1,'first');
            obs_percentiles(gg,g) = gsizel(tt);
            [~,obs_crit_perc(gg,g)] = pmsoulsby(obs_percentiles(gg,g),1);
        end

        
        for mm = 1:length(mlists)
            tSeries2 = R(gg)*tSeries(ismember(gTime(:,2),mlists{mm}));
        obs_model_percentiles(gg,mm,g) = 100*...
            (length(find(tSeries2 >= obs_crit_perc(gg,g)))./...
            length(find(~isnan(tSeries2))));
        obs_perc_over_pt1(gg,mm) = 100* ...
            (length(find(tSeries2 >= 0.1))) ./...
            length(tSeries2);
                obs_perc_over_pt2(gg,mm) = 100* ...
            (length(find(tSeries2 >= 0.2))) ./...
            length(tSeries2);
        end
    end
    clear gall Tseries tSeries tSeries2
end
clear nc

clear ii jj gg g var tt mc mp gsizel box_obs_in ans
obs_crit_perc = squeeze(obs_crit_perc);
obs_model_percentiles = squeeze(obs_model_percentiles);
obs_percentiles = squeeze(obs_percentiles);
clear frac
save compare_sediment_percentiles_to_skinstress_tseries_cohesive_waves