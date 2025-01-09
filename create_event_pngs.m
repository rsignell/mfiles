%create_event_pngs

ncclear
mp.dir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';
mp.file = 'analysis2_events.mat';

%Creates the plots for the recurrence interval stuff

%% Bed mobility
% mp.pvar = 'perc_over(1,1,:,:)'; %<------CHANGE THIS (thresh,timeperiod,lon,lat) 
% mp.title = '% of Time Above Threshold, 0.1 Pa'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 100];
% mp.tpoints = [0:20:100];
% mp.print = 'mab_yrly_pover_pt1.png'; %<--------CHANGE THIS

%% Event analysis
% mp.pvar = 'med_event_dur(1,1,:,:)'; %<------CHANGE THIS (thresh,timeperiod,lon,lat) 
% mp.title = 'Median Event Duration (Days); Threshold = 0.1 Pa'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 2];
% mp.tpoints = [0:0.5:2];
% mp.print = 'mab_yrly_meded_pt1.png'; %<--------CHANGE THIS

%% Quiescent period analysis
mp.pvar = 'mean_quis_dur(1,3,:,:)'; %<------CHANGE THIS (thresh,timeperiod,lon,lat) 
mp.title = 'Mean Quiescent Duration, Spring (Days); Threshold = 0.1 Pa'; %<--------CHANGE THIS
mp.plottype = 'linear';    %Log or linear
mp.plotclean = 'esp_plots';
mp.mrange = [0 1];
mp.tpoints = [0:0.25:1];
mp.print = 'mab_mean_spr_qu_pt1.png'; %<--------CHANGE THIS

%% mid range values analysis
% mp.pvar = 'ptween(1,3,:,:)'; %<------CHANGE THIS (thresh,timeperiod,lon,lat) 
% mp.title = '% of Time 0.1 Pa <= \tau_w_c <= 0.5 Pa, Spring'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 100];
% mp.tpoints = [0:20:100];
% mp.print = 'mab_spring_ptween.png'; %<--------CHANGE THIS
% 
% mp.pvar = 'ptween(2,1,:,:)'; %<------CHANGE THIS (thresh,timeperiod,lon,lat) 
% mp.title = '% of Time 0.5 Pa <= \tau_w_c <= 1 Pa'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 20];
% mp.tpoints = [0:5:20];
% mp.print = 'mab_big_ptween.png'; %<--------CHANGE THIS

%% Recurr Interval
% mp.file = 'analysis2_events.mat';
% mp.pvar = 'mean_rec_interval(1,1,:,:)'; %<------CHANGE THIS (thresh,timeperiod,lon,lat) 
% mp.title = 'Recurrence Interval (Days), 0.1 Pa'; %<--------CHANGE THIS
% mp.plottype = 'linear_flip';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 10];
% mp.tpoints = [0:2:10];
% mp.print = 'mab_yrly_ri_pt1.png'; %<--------CHANGE THIS





%% Load up everything
%Data
cd(mp.dir)
load(fullfile(mp.dir,mp.file))

%Bathymetry
nc1 = mDataset('http://geoport.whoi.edu/thredds/dodsC/bathy/etopo2_v2c.nc');
lon1 = nc1{'lon'}(:);
lat1 = nc1{'lat'}(:);
jj = find((lon1 >= min(min(lon))) & (lon1 <= max(max(lon))));
ii = find((lat1 >= min(min(lat))) & (lat1 <= max(max(lat))));
lat1 = lat1(ii); lon1 = lon1(jj);
topo = nc1{'topo'}(ii,jj);
close(nc1), clear nc1 ii jj

%Coastline, states
PBFile = fullfile('C:\Users\sdalyander\Documents\Maps\WDB','namer-pby.txt');
PB = read_wdbfile(PBFile);
PBFile = fullfile('C:\Users\sdalyander\Documents\Maps\WDB','namer-cil.txt');
PB2 = read_wdbfile(PBFile);


%% Analysis
pt = squeeze(perc_over(1,1,:,:));
sprlow = squeeze(mean_quis_dur(1,3,:,:));
good = pt >= 50 & sprlow >= 0.75;
%% Plots
%Feed in data
eval(['myData = squeeze(' mp.pvar ');'])
myData(myData == Inf) = NaN;    %Recurrence interval
close all
figure
clear mAx2
if strcmp(mp.plottype,'log')
    [mAx,mAx2,cbar] = logpcolorpsd(lon, lat, myData,mp.mrange,0);
elseif strcmp(mp.plottype,'linear')
    [mAx,cbar] = pcolorpsd(lon,lat,myData,mp.mrange,mp.tpoints);
elseif strcmp(mp.plottype,'linear_flip')
    [mAx,cbar] = pcolorpsd_flip(lon,lat,myData,mp.mrange,mp.tpoints);
elseif strcmp(mp.plottype,'linear_flip_cut')
    [mAx,cbar] = pcolorpsd_flip_cutoff(lon,lat,myData,mp.mrange,mp.tpoints);
elseif strcmp(mp.plottype,'linear_cut')
    [mAx,cbar] = pcolorpsd_cutoff(lon,lat,myData,mp.mrange,mp.tpoints);
end
cBarLabel = mp.title;
esp_plots
 print('-dpng',fullfile(mp.dir,mp.print))
%  close
