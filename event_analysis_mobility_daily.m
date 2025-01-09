% bed_mobility_by_sed_core_percentage

ncclear, clc
load perc_grainsize_crit_smaller_percentiles box_obs_in mc gsizel p

split = find(gsizel < 0.004e-3,1,'last');

gList = 10:10:90;

%analysis time periods
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

perc_days_events = NaN([size(mc,1) length(mlists) size(gList)]);
avg_days_between_events = NaN([size(mc,1) length(mlists) size(gList)]);
obs_crit_perc = NaN([size(mc,1) size(gList)]);
obs_percentiles = NaN([size(mc,1) size(gList)]);

mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
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
    total = total*R(gg);
    
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
        
        isOver = total >= obs_crit_perc(gg,g);
        %Pad to an even day
        isOver = [isOver(:)' zeros(1,24-rem(length(isOver),24))];
        newTime = [time(:)' zeros(1,24-rem(length(time),24))];
        isOver = reshape(isOver,24,numel(isOver)/24);
        newTime = reshape(newTime,24,numel(newTime)/24);
        dayList = newTime(1,:);
        gTimeDay = gregorian(dayList);
        isOver = sum(isOver,1);
        isOver = isOver > 0;
        
        for mm = 1:length(mlists)
            inThis = ismember(gTimeDay(:,2),mlists{mm});
            perc_days_events(gg,mm,g) = sum(isOver(inThis))/length(isOver(inThis));
            
            %Quiscent periods
            dQuis = diff([0 ~isOver(inThis) 0]);
            qStart = find(dQuis == 1);
            qEnd = find(dQuis == -1)-1;
            qlength = qEnd - qStart;
            
            avg_days_between_events(gg,mm,g) = nanmean(qlength);
        end
        clear sStart qStart slength qlength
        
    end
    clear gall total tide wave curr res
end
clear nc

clear ii jj gg g tt mc mp gsizel box_obs_in ans clear tauwc tauc_tide
obs_crit_perc = squeeze(obs_crit_perc);
clear tauc tauc_res tauw split highfrac lowfrac low high f hmss newStress pow

save mobility_event_analysis_daily_reducestress10percent