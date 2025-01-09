% bed_mobility_by_sed_core_percentage

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

mean_event_dur = NaN([NL length(mlists)]);
mean_quis_dur = NaN([NL length(mlists)]);
std_event_dur = NaN([NL length(mlists)]);
std_quis_dur = NaN([NL length(mlists)]);
mean_rec_interval = NaN([NL length(mlists)]);


cuts = 1:500:NL-500;
cuts2 = 500:500:NL;

for cc = 1:length(cuts)
    tSeries = nc{'tauwc'}(:,cuts(cc):cuts2(cc));
    crit_stress = nc{'tau_crit'}(cuts(cc):cuts2(cc));
    
    if all(isnan(crit_stress)), continue, end
    disp(['On cc = ' num2str(cc) ' of ' num2str(length(cuts)) '.'])
      
    for gg = 1:length(crit_stress)
        if isnan(crit_stress(gg)),continue,end
        total = tSeries(:,gg);
        
        dStorm = diff([0; total >= crit_stress(gg); 0]);
        sStart = time(dStorm == 1);
        sEnd = time(find(dStorm == -1)-1);
        slength = sEnd - sStart;
        sStart(slength <= (1.5/24)) = [];
        slength(slength <= (1.5/25)) = []; %Single point over
        clear sEnd
        
        %Quiscent periods
        dQuis = diff([0; total < crit_stress(gg); 0]);
        qStart = time(dQuis == 1);
        qEnd = time(find(dQuis == -1)-1);
        qlength = qEnd - qStart;
        qStart(qlength <= (1.5/24)) = [];
        qlength(qlength <= (1.5/24)) = []; %Single point over
        clear qEnd
        
        sStart = datevec(sStart);
        qStart = datevec(qStart);
        
        for mm = 1:length(mlists)
            
            mean_event_dur(cuts(cc)+gg-1,mm) = ...
                nanmean(slength(ismember(sStart(:,2),mlists{mm})));
            mean_quis_dur(cuts(cc)+gg-1,mm) = ...
                nanmean(qlength(ismember(qStart(:,2),mlists{mm})));
            mean_rec_interval(cuts(cc)+gg-1,mm) = ...
                (1/24)*(length(find(ismember(gTime(:,2),mlists{mm})))./...
                length(find(ismember(sStart(:,2),mlists{mm}))));
            
        end
        clear sStart qStart slength qlength
        
    end
    clear gall total tide wave curr res
end
clear nc

clear ii jj gg g tt mc mp gsizel box_obs_in ans clear tauwc tauc_tide
clear tauc tauc_res tauw split highfrac lowfrac low high f hmss newStress pow

save mobility_event_analysis