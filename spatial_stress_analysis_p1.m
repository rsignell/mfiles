%% spatial_stress_analysis_p1

%Includes general statistical analysis (mean, std by stress type and
%season) and geologic basic analysis (bed mobility for phi class, fraction
%of wave/current/tide induced stress when bed is mobil)
ncclear
tic
%% Where data is from, where everything will be saved to
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SABGOM_SWAN7/sabgom_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\SABGOM';

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

%Set 1:  basic stats
set1 = {'gmean'; 'gstd'; 'mmean'; 'mstd'; 'mmed'; 'm95'; 'm05'; 'm75'; 'm25'; 'm10'; 'm20'; ...
    'm30'; 'm40'; 'm60'; 'm70'; 'm80'; 'm90';'m16';'m84'};
com1s = {'nangeomean('; 'nangeostd('; 'nanmean('; 'nanstd('; 'nanmedian('; 'prctile('; 'prctile('; ...
    'prctile('; 'prctile('; 'prctile('; 'prctile('; 'prctile('; 'prctile('; ...
    'prctile('; 'prctile('; 'prctile('; 'prctile('; 'prctile('; 'prctile('; 'prctile('};
com1e = {');';');';');'; ');'; ');'; ',95);'; ',5);'; ',75);'; ',25);'; ',10);'; ',20);'; ...
    ',30);'; ',40);'; ',60);'; ',70);'; ',80);'; ',90);'; ',16);'; ',84);'};

%Set 2:  geologic thresholds
gsize = [62e-6 1e-3 2e-3];


%% Set up empties
%Set 1...[tau_type(x5) timeframe(x5) spatial]
gmean = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
gstd = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
mmean = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
mstd = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
mmed = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m95 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m84  = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m05 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m75 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m25 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m10 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m16 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m20 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m30 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m40 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m60 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m70 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m80 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
m90 = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);

%Set 2...[gsize timeframe(x5) spatial]
mbed_mobil = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);
tauw_perc = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);
tauc_perc = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);
taut_perc = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);
taur_perc = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);

taut_percc = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);
taur_percc = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);

taut2s = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);


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
        tauw = nc.data('tauw',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        tauc = nc.data('tauc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        taut = nc.data('tauc_tide',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        taur = nc.data('tauc_res',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        disp('Analyzing data')
        
        tauwc(tauwc==0) = NaN;
        
        tauw_wc = tauw./tauwc;
        tauc_wc = tauc./tauwc;
        taut_wc = taut./tauwc;
        taur_wc = taur./tauwc;
        
        tauc(tauc == 0) = NaN;
        taut_c = taut./tauc;
        taur_c = taur./tauc;
        
        tau_s = tauw + taur;
        tau_s(tau_s == 0) = NaN;
        ratiot_s = taut./tau_s;
        
        %Cycle through
        for mm = 1:length(mlists)
            inThis = find(ismember(gtime(:,2),mlists{mm}));
            
            for tt = 1:length(ttypes)   %w, c, wc
                %Set 1:  General stats
                for ss = 1:length(set1)
                    eval([set1{ss} '(tt,mm,jt,it)=' com1s{ss} 'tau' ttypes{tt} ...
                        '(inThis,:,:)' com1e{ss}])
                end
            end
            
            %Set 2:  Geologic
            for gg = 1:length(gsize)
                [mbed_mobil(gg,mm,jt,it),tau_crit] = ...
                    pmsoulsby(gsize(gg),tauwc(inThis,:,:));
                overCheck = double((tauwc >= tau_crit));
                overCheck(overCheck == 0) = NaN;
                tauw_perc(gg,mm,jt,it) = nanmean(overCheck(inThis,:,:).*tauw_wc(inThis,:,:),1);
                tauc_perc(gg,mm,jt,it) = nanmean(overCheck(inThis,:,:).*tauc_wc(inThis,:,:),1);
                taut_perc(gg,mm,jt,it) = nanmean(overCheck(inThis,:,:).*taut_wc(inThis,:,:),1);
                taur_perc(gg,mm,jt,it) = nanmean(overCheck(inThis,:,:).*taur_wc(inThis,:,:),1);
                
                taut_percc(gg,mm,jt,it) = nanmean(overCheck(inThis,:,:).*taut_c(inThis,:,:),1);
                taur_percc(gg,mm,jt,it) = nanmean(overCheck(inThis,:,:).*taur_c(inThis,:,:),1);
                
                taut2s(gg,mm,jt,it) = nanmean(overCheck(inThis,:,:).*ratiot_s(inThis,:,:),1);
            end
            
            tauw_perc(gg+1,mm,jt,it) = nanmean(tauw_wc(inThis,:,:),1);
            tauc_perc(gg+1,mm,jt,it) = nanmean(tauc_wc(inThis,:,:),1);
            taut_perc(gg+1,mm,jt,it) = nanmean(taut_wc(inThis,:,:),1);
            taur_perc(gg+1,mm,jt,it) = nanmean(taur_wc(inThis,:,:),1);
            
            taut_percc(gg+1,mm,jt,it) = nanmean(taut_c(inThis,:,:),1);
            taur_percc(gg+1,mm,jt,it) = nanmean(taur_c(inThis,:,:),1);
            
            taut2s(gg+1,mm,jt,it) = nanmean(ratiot_s(inThis,:,:),1);
        end
    end
    clear tauwc tauw tauc tauw_wc tauc_wc taut_wc taur_wc taut taur
end

%Spread
mspread = m95 - m05;
% close(nc), clear nc
clear nc
clear overCheck com1e com1s down_dur event_dur set1 ss tau_crit test tt
clear gg gtime ii iiL inThis it jj jjL jt mm mp
disp('Be sure to save results!!')
toc