%% spatial_stress_analysis_p7

%Low-passing the outcome

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

%stress types, analysis time periods
ttypes = {'wc'; 'w'; 'c'; 't'; 'r'};
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};


%% Set up empties [tau timeframe(x5) spatial]
lp_std = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);    %mean event duration
hf_std = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
mean_all = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);

%Assumes hourly data!!
%% Run the loop
for jj = 1:length(jjL)-1
    for ii = 1:length(iiL)-1
        disp(['Data set ' num2str(jj) ',' num2str(ii) ' of ' ...
            num2str(length(jjL)-1) ',' num2str(length(iiL)-1)])
        
        jt = jjL(jj):jjL(jj+1)-1;
        it = iiL(ii):iiL(ii+1)-1;
        
        testCut = test(jt,it);
        
        disp('Loading data')
        %         tauwc = squeeze(double(nc{'tauwc'}(:,jt,it)));
        %         tauw = squeeze(double(nc{'tauw'}(:,jt,it)));
        %         tauc = squeeze(double(nc{'tauc'}(:,jt,it)));
                        clear tauwc tauw tauc tauw_wc tauc_wc taut_wc taur_wc taut taur
        tauwc = nc.data('tauwc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        tauw = nc.data('tauw',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        tauc = nc.data('tauc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        taut = nc.data('tauc_tide',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        taur = nc.data('tauc_res',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        disp('Analyzing data')
        
        for jj2 = 1:length(jt)
            for ii2 = 1:length(it)
                if isnan(testCut(jj2,ii2)), continue, end
                for tt = 1:length(ttypes)
                    
                    
                    %Cycle through
                    for mm = 1:length(mlists)
                        inThis = find(ismember(gtime(:,2),mlists{mm}));
                        eval(['tau_this = tau' ttypes{tt} '(inThis,jj2,ii2);'])
                        mean_all(tt,mm,jt(jj2),it(ii2)) = nanmean(tau_this);
                        tau_this = tau_this - nanmean(tau_this);
                        [tau_lp,jdlp] = plfilt(tau_this,time(inThis),1);
                        tau_lp = smart_interp(jdlp,tau_lp,time(inThis),7);
                        tau_hp = tau_this - tau_lp;
                        lp_std(tt,mm,jt(jj2),it(ii2)) = nanstd(tau_lp);
                        hf_std(tt,mm,jt(jj2),it(ii2)) = nanstd(tau_hp);
                        
                        
                    end
                end
            end
        end
    end
end
% close(nc), clear nc
clear nc tau_lp tau_hp tau_this jdlp inThis
clear overCheck com1e com1s down_dur event_dur set1 ss tau_crit test tt cc
clear gg gtime ii iiL inThis it jj jjL jt mm testCut dQuis dStorm ii2 jj2 mp ttau
disp('Be sure to save results!!')
toc