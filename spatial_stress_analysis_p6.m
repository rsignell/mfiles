%% spatial_stress_analysis_p2

%Event-based analysis

%IMPORTANTS NOTES:  Time series needs to be continuous, and for the
%recurrence interval analysis assumes total time span of < 1 year
%Storms are put into the season in which the start time falls

ncclear, clc

tic
%% Where data is from, where everything will be saved to
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';


%% Open the dataset
% nc = mDataset(mp.urlH);
% time = nj_time(nc,'tauwc');
% lon = double(nc{'lon'}(:));
% lat = double(nc{'lat'}(:));
nc = ncdataset(mp.urlH);
time = nc.time('time');
lon = nc.data('lon');
lat = nc.data('lat');

gtime = datevec(time);

%Pull down one time step as a mask
% test = nc{'tauwc'}(100,:,:);
test = squeeze(double(nc.data('tauwc',[100 1 1],[100 size(lon,1) size(lon,2)])));

%Set up divisions
nn = 50;
jjL = unique([1:nn:size(test,1) size(test,1)]);
iiL = unique([1:nn:size(test,2) size(test,2)]);
clear nn

%analysis time periods
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

%Set 3:  Event based

load perc_grainsize_crit_smaller_percentiles box_obs_in cstress_grains

%% Set up empties [cvalue timeframe(x5) spatial]
mean_event_dur = NaN([size(cstress_grains,2) length(mlists) size(lon,1) size(lon,2)]);    %mean event duration
mean_quis_dur = NaN([size(cstress_grains,2) length(mlists) size(lon,1) size(lon,2)]); %mean duration between events
mean_rec_interval = NaN([size(cstress_grains,2) length(mlists) size(lon,1) size(lon,2)]); %time period/events, time period in days
med_event_dur = NaN([size(cstress_grains,2) length(mlists) size(lon,1) size(lon,2)]);
med_quis_dur = NaN([size(cstress_grains,2) length(mlists) size(lon,1) size(lon,2)]);
perc_over = NaN([size(cstress_grains,2) length(mlists) size(lon,1) size(lon,2)]);
ptween = NaN([size(cstress_grains,2)-1 length(mlists) size(lon,1) size(lon,2)]);


%Assumes hourly data!!
%% Run the loop
for jj = 1:length(jjL)-1
    for ii = 1:length(iiL)-1
        disp(['Data set ' num2str(jj) ',' num2str(ii) ' of ' ...
            num2str(length(jjL)-1) ',' num2str(length(iiL)-1)])
        
        jt = jjL(jj):jjL(jj+1)-1;
        it = iiL(ii):iiL(ii+1)-1;
        
        disp('Loading data')
        %         tauwc = squeeze(double(nc{'tauwc'}(:,jt,it)));
        %         tauw = squeeze(double(nc{'tauw'}(:,jt,it)));
        %         tauc = squeeze(double(nc{'tauc'}(:,jt,it)));
        tauwc = nc.data('tauwc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        disp('Analyzing data')
        
        testCut = test(jt,it);
        
        %Cycle through
        %Point by point analysis
        for jj2 = 1:length(jt)
            for ii2 = 1:length(it)
                if isnan(testCut(jj2,ii2)), continue, end
                
                inThis = find(box_obs_in(:,1) == jt(jj2) & box_obs_in(:,2) == it(ii2));
                if isempty(inThis), continue, end
                ttau = tauwc(:,jj2,ii2);
                for gg = 1:size(cstress_grains,2)                   
                    
                    critvalue = nanmean(cstress_grains(inThis,gg));
                    
                    %Storms
                    dStorm = diff([0; ttau >= critvalue; 0]);
                    sStart = time(dStorm == 1);
                    sEnd = time(find(dStorm == -1)-1);
                    slength = sEnd - sStart;
                    sStart(slength <= (1.5/24)) = [];
                    slength(slength <= (1.5/25)) = []; %Single point over
                    clear sEnd
                    
                    %Quiscent periods
                    dQuis = diff([0; ttau < critvalue; 0]);
                    qStart = time(dQuis == 1);
                    qEnd = time(find(dQuis == -1)-1);
                    qlength = qEnd - qStart;
                    qStart(qlength <= (1.5/24)) = [];
                    qlength(qlength <= (1.5/24)) = []; %Single point over
                    clear qEnd
                    
                    sStart = datevec(sStart);
                    qStart = datevec(qStart);
                    
                    for mm = 1:length(mlists)
                        if gg < size(cstress_grains,2)
                            critplus1 = nanmean(cstress_grains(inThis,gg+1));
                            ptween(gg,mm,jt(jj2),it(ii2)) = ...
                                100*(length(find(ttau >= critvalue & ...
                                ttau < critplus1 & ...
                                ismember(gtime(:,2),mlists{mm})))./ ...
                                length(find(ismember(gtime(:,2),mlists{mm}))));
                        end
                        
                        mean_event_dur(gg,mm,jt(jj2),it(ii2)) = ...
                            nanmean(slength(ismember(sStart(:,2),mlists{mm})));
                        mean_quis_dur(gg,mm,jt(jj2),it(ii2)) = ...
                            nanmean(qlength(ismember(qStart(:,2),mlists{mm})));
                        mean_rec_interval(gg,mm,jt(jj2),it(ii2)) = ...
                            (1/24)*(length(find(ismember(gtime(:,2),mlists{mm})))./...
                            length(find(ismember(sStart(:,2),mlists{mm}))));
                        med_event_dur(gg,mm,jt(jj2),it(ii2)) = ...
                            nanmedian(slength(ismember(sStart(:,2),mlists{mm})));
                        med_quis_dur(gg,mm,jt(jj2),it(ii2)) = ...
                            nanmedian(qlength(ismember(qStart(:,2),mlists{mm})));
                        
                        perc_over(gg,mm,jt(jj2),it(ii2)) = ...
                            100*(length(find(ttau >= critvalue & ...
                            ismember(gtime(:,2),mlists{mm})))./ ...
                            length(find(ismember(gtime(:,2),mlists{mm}))));
                    end
                    clear sStart qStart slength qlength
                end
            end
        end
    end
    clear tauwc tauw tauc tauw_wc tauc_wc
end
% close(nc), clear nc
clear nc
clear overCheck com1e com1s down_dur event_dur set1 ss tau_crit test tt cc
clear gg gtime ii iiL inThis it jj jjL jt mm testCut dQuis dStorm ii2 jj2 mp ttau
clear bads ans critplus1 critvalue first fover isIn last m n sum_mc x y
disp('Be sure to save results!!')
toc