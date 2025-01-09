%% spatial_stress_analysis_p1

%Includes general statistical analysis (mean, std by stress type and
%season) and geologic basic analysis (bed mobility for phi class, fraction
%of wave/current/tide induced stress when bed is mobil)
ncclear
tic
%% Where data is from, where everything will be saved to
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CompGrid/swan_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\WavesOnly';

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
test = squeeze(double(nc.data('tauw',[100 1 1],[100 size(lon,1) size(lon,2)])));

%Set up divisions
nn = 50;
jjL = unique([1:nn:size(test,1) size(test,1)]);
iiL = unique([1:nn:size(test,2) size(test,2)]);
clear nn


%stress types, analysis time periods
% ttypes = {'wc'; 'w'; 'c'; 't'; 'r'}; %Waves only
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

%Set 1:  basic stats
set1 = {'mmean'; 'mstd'; 'mmed'; 'm95'; 'm05'; 'm75'; 'm25'; 'm10'; 'm20'; ...
    'm30'; 'm40'; 'm60'; 'm70'; 'm80'; 'm90'};
com1s = {'nanmean('; 'nanstd('; 'nanmedian('; 'prctile('; 'prctile('; ...
    'prctile('; 'prctile('; 'prctile('; 'prctile('; 'prctile('; 'prctile('; ...
    'prctile('; 'prctile('; 'prctile('; 'prctile('; 'prctile('};
com1e = {');'; ');'; ');'; ',95);'; ',5);'; ',75);'; ',25);'; ',10);'; ',20);'; ...
    ',30);'; ',40);'; ',60);'; ',70);'; ',80);'; ',90);'};

%Set 2:  geologic thresholds
gsize = [62e-6 1e-3 2e-3];


%% Set up empties
%Set 1...[tau_type(x5) timeframe(x5) spatial]
mmean = NaN([length(mlists) size(lon,1) size(lon,2)]);
mstd = NaN([length(mlists) size(lon,1) size(lon,2)]);
mmed = NaN([length(mlists) size(lon,1) size(lon,2)]);
m95 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m05 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m75 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m25 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m10 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m20 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m30 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m40 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m60 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m70 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m80 = NaN([length(mlists) size(lon,1) size(lon,2)]);
m90 = NaN([length(mlists) size(lon,1) size(lon,2)]);

%Set 2...[gsize timeframe(x5) spatial]
mbed_mobil = NaN([length(gsize) length(mlists) size(lon,1) size(lon,2)]);

%% Run the loop
for jj = 1:length(jjL)-1
    for ii = 1:length(iiL)-1
        disp(['Data set ' num2str(jj) ',' num2str(ii) ' of ' ...
            num2str(length(jjL)-1) ',' num2str(length(iiL)-1)])
        
        jt = jjL(jj):jjL(jj+1)-1;
        it = iiL(ii):iiL(ii+1)-1;
        
        disp('Loading data')
        
        tauw = nc.data('tauw',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        
        disp('Analyzing data')
        
        
        %Cycle through
        for mm = 1:length(mlists)
            
            %Set 1:  General stats
            for ss = 1:length(set1)
                eval([set1{ss} '(mm,jt,it)=' com1s{ss} 'tauw' ...
                    '(ismember(gtime(:,2),mlists{mm}),:,:)' com1e{ss}])
            end
            
            %Set 2:  Geologic
            for gg = 1:length(gsize)
                [mbed_mobil(gg,mm,jt,it)] = ...
                    pmsoulsby(gsize(gg),tauw(ismember(gtime(:,2),mlists{mm}),:,:));

            end
            
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