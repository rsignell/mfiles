%% spatial_stress_analysis

clc
ncclear
%% Where data is from, where everything will be saved to
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';

%% Open the dataset
nc = mDataset(mp.urlH);
time = nj_time(nc,'tauwc');
lon = double(nc{'lon'}(:));
lat = double(nc{'lat'}(:));
% nc = ncdataset(mp.urlH);
% time = nc.time('time');
% lon = nc.data('lon');
% lat = nc.data('lat');

gtime = datevec(time);

%Pull down one time step as a mask
test = nc{'tauwc'}(100,:,:);

% test = squeeze(double(nc.data('tauwc',[100 1 1],[100 size(lon,1) size(lon,2)])));
[jjL,iiL] = find(~isnan(test));
clear test

%stress types, analysis time periods
ttypes = {'w'; 'c'; 'wc'};
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

%Set 1:  basic stats
set1 = {'mmean'; 'mstd'; 'm95'};
com1s = {'nanmean('; 'nanstd('; 'prctile('};
com1e = {');'; ');'; ',95);'};

%Set 2:  geologic thresholds
gsize = [62e-6 2e-3];

%Set 3:  event based (biological?)...not on yet
% cvalues = [0.1 0.5 1];  %Stress values


%% Set up empties
%Set 1...[tau_type(x3) timeframe(x5) spatial]
mmean = NaN([3 length(mlists) size(lon,1) size(lon,2)]);
mstd = NaN([3 length(mlists) size(lon,1) size(lon,2)]);
m95 = NaN([3 length(mlists) size(lon,1) size(lon,2)]);

%Set 2...[gsize timeframe(x5) spatial]
mbed_mobil = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);
tauw_perc = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);
tauc_perc = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);

% %Set 3....[cvalue timeframe(x5) spatial]
% event_dur = NaN([length(cvalues) length(mlists) size(lon,1) size(lon,2)]);    %mean event duration
% down_dur = NaN([length(cvalues) length(mlists) size(lon,1) size(lon,2)]); %mean duration between events

%% Run the loop
for ind = 1:length(jjL)
    disp(['Data point ' num2str(ind) ' of ' num2str(length(jjL))])
    
    jj=jjL(ind);
    ii=iiL(ind);
    
%     %Kludgy bit...crashes without it
%     if rem(ind,10)==0
%         disp('Pausing')
%         pause(2*60)
%     end
    
    tauwc = squeeze(double(nc{'tauwc'}(:,jj,ii)));
    tauw = squeeze(double(nc{'tauw'}(:,jj,ii)));
    tauc = squeeze(double(nc{'tauc'}(:,jj,ii)));
    
    %
    %     tauwc = nc.data('tauwc',[1 jj ii],[length(time) jj ii]);
    %     tauw = nc.data('tauw',[1 jj ii],[length(time) jj ii]);
    %     tauc = nc.data('tauc',[1 jj ii],[length(time) jj ii]);
    
    tauw_wc = tauw./tauwc;
    tauc_wc = tauc./tauwc;
    
    %Cycle through
    for mm = 1:length(mlists)
        inThis = find(ismember(gtime(:,2),mlists{mm}));
        
        for tt = 1%:length(ttypes)   %w, c, wc
            %Set 1:  General stats
            for ss = 1%:length(set1)
                eval([set1{ss} '(tt,mm,jj,ii)=' com1s{ss} 'tau' ttypes{tt} ...
                    '(inThis)' com1e{ss}])
            end
        end
        
%         %Set 2:  Geologic
%         for gg = 1:length(gsize)
%             [mbed_mobil(gg,mm,jj,ii),tau_crit] = pmsoulsby(gsize(gg),tauwc(inThis));
%             inThisO = find(ismember(gtime(:,2),mlists{mm}) & (tauwc >= tau_crit));
%             tauw_perc(gg,mm,jj,ii) = nanmean(tauw_wc(inThisO));
%             tauc_perc(gg,mm,jj,ii) = nanmean(tauc_wc(inThisO));
%         end
    end
    
    clear tauwc tauw tauc tauw_wc tauc_wc
end
close(nc)
%clear nc
