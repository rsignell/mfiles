% bed_mobility_by_sed_core_percentage

ncclear, clc
load compare_sediment_percentiles_to_skinstress_tseries_cohesive R box_obs_in mc gsizel p

split = find(gsizel < 0.004e-3,1,'last');

gList = 10:10:90;

%analysis time periods
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

mean_event_dur = NaN([size(mc,1) length(mlists) size(gList)]);
mean_quis_dur = NaN([size(mc,1) length(mlists) size(gList)]);
std_event_dur = NaN([size(mc,1) length(mlists) size(gList)]);
std_quis_dur = NaN([size(mc,1) length(mlists) size(gList)]);
mean_rec_interval = NaN([size(mc,1) length(mlists) size(gList)]);
obs_crit_perc = NaN([size(mc,1) size(gList)]);
obs_percentiles = NaN([size(mc,1) size(gList)]);

mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau_zm5_bdir/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';
nc = cfdataset(mp.urlH);
time = nc.time('time');
gtime = datevec(time);
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
    %skin friction
    total = R(gg)*total;
    
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
        
        dStorm = diff([0; total >= obs_crit_perc(gg,g); 0]);
        sStart = time(find(dStorm == 1));
        sEnd = time(find(dStorm == -1)-1);
        slength = sEnd - sStart;
        sStart(slength <= (1.5/24)) = [];
        slength(slength <= (1.5/25)) = []; %Single point over
        clear sEnd
        
        %Quiscent periods
        dQuis = diff([0; total < obs_crit_perc(gg,g); 0]);
        qStart = time(find(dQuis == 1));
        qEnd = time(find(dQuis == -1)-1);
        qlength = qEnd - qStart;
        qStart(qlength <= (1.5/24)) = [];
        qlength(qlength <= (1.5/24)) = []; %Single point over
        clear qEnd
        
        sStart = datevec(sStart);
        qStart = datevec(qStart);
        
        for mm = 1:length(mlists)
            
            mean_event_dur(gg,mm,g) = ...
                nanmean(slength(ismember(sStart(:,2),mlists{mm})));
            mean_quis_dur(gg,mm,g) = ...
                nanmean(qlength(ismember(qStart(:,2),mlists{mm})));
            mean_rec_interval(gg,mm,g) = ...
                (1/24)*(length(find(ismember(gtime(:,2),mlists{mm})))./...
                length(find(ismember(sStart(:,2),mlists{mm}))));
            
        end
        clear sStart qStart slength qlength
        
    end
    clear gall total tide wave curr res
end
clear nc

clear ii jj gg g tt mc mp gsizel box_obs_in ans clear tauwc tauc_tide
obs_crit_perc = squeeze(obs_crit_perc);
clear tauc tauc_res tauw split highfrac lowfrac low high f hmss newStress pow

save mobility_event_analysis