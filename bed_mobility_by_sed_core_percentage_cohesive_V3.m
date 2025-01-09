% bed_mobility_by_sed_core_percentage
%%
ncclear, clc

mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SABGOM_SWAN7_point/sabgom_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\SABGOM\stress_analysis';


mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

nc = mDataset(mp.urlH);
time = nj_time(nc,'tauwc');
twc = nc{'tauwc'}(1,:);
gTime = datevec(time);
NL = length(twc);
locs = [nc{'lon'}(:) nc{'lat'}(:)];

obs_model_percentiles = NaN(NL,length(mlists));

cuts = 1:500:NL-500;
cuts2 = 500:500:NL;

%%
for cc = 1:length(cuts)
    tSeries = nc{'tauwc'}(:,cuts(cc):cuts2(cc));
    crit_stress = nc{'tau_crit'}(cuts(cc):cuts2(cc));
    
    if all(isnan(crit_stress)), continue, end
    disp(['On cc = ' num2str(cc) ' of ' num2str(length(cuts)) '.'])
    
%     tSeries = nc{'tauwc'}(:,gg);
%     crit_stress = nc{'tau_crit'}(gg);
    crit_stress = repmat(crit_stress(:)',[size(tSeries,1) 1]);
    isOver = tSeries >= crit_stress;
    isNaN = ~isnan(tSeries);
        for mm = 1:length(mlists)
            isIn = ismember(gTime(:,2),mlists{mm});
            checkMat = 100.*sum(isOver(isIn,:),1)./sum(isNaN(isIn,:),1);
            obs_model_percentiles(cuts(cc):cuts2(cc),mm) = checkMat;
        end
    clear tSeries tSeries2
end
crit_stress = nc{'tau_crit'}(:);
obs_model_percentiles(isnan(crit_stress),:) = NaN;

%%
clear nc

% clear isNaN isOver isIn cuts cuts2 crit_stress checkMat cc NL num time twc gTime mm
% clear ii jj gg g var tt mc mp gsizel box_obs_in ans
save compare_sediment_percentiles_to_skinstress_tseries_cohesive