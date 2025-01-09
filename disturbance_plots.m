%% disturbance_plots

ncclear,clc
%%
mask = load('analysis1_gen_geo', 'mmean');
mask = squeeze(mask.mmean(1,1,:,:));
mask(mask >=0) = 1;

mp.dir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';
cd(mp.dir)

%% ("Actual") Bed mobility
% mp.file = 'analysis4_disturbance.mat';
% mp.pvar = 'mbed_mobil_act(1,1,:,:)'; %<------CHANGE THIS (stress mechanism,timeperiod,lon,lat) 
% mp.title = 'Yearly Bed Mobility, Actual Grain Size'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 100];
% mp.tpoints = [0:25:100];
% mp.print = 'mab_yrly_bm_act.png'; %<--------CHANGE THIS

%% Percent below threshold
mp.file = 'perc_grainsize_crit_smaller_percentiles.mat';
mp.pvar = 'perc_crit_smaller_obs(5,:,:)'; %<------CHANGE THIS (percentile,lon,lat) 
mp.title = {'Percentage of Sediment'; '\tau_c_r_i_t < 50th Percentile Stress'}; %<--------CHANGE THIS
mp.plottype = 'linear';    %Log or linear
mp.plotclean = 'esp_plots';
mp.mrange = [0 100];
mp.tpoints = [0:25:100];
mp.print = 'mab_perc_grain_crit_below_p50.png'; %<--------CHANGE THIS

%% event driven
% mp.file = 'analysis6_events_obs_sediment.mat';
% mp.pvar = 'mean_quis_dur(5,1,:,:)'; %<------CHANGE THIS (sed. perc,timeperiod,lon,lat) 
% mp.title = {'Yearly Average Quiscient Duration'; '50th Percentile of Sediment'}; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 10];
% mp.tpoints = [0:2:10];
% mp.print = 'mab_yrly_quis_act_50p.png'; %<--------CHANGE THIS



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


%% Plots
%Feed in data
eval(['myData = squeeze(' mp.pvar ');'])
myData(myData == Inf) = NaN;    %Recurrence interval
close all
figure, orient landscape
subplot(1,2,1)
clear mAx2
if strcmp(mp.plottype,'log')
    [mAx,mAx2,cbar] = logpcolorpsd(lon, lat, myData,mp.mrange,0);
elseif strcmp(mp.plottype,'linear')
    [mAx,cbar] = pcolorpsd(lon,lat,myData,mp.mrange,mp.tpoints);
end
cBarLabel = mp.title;
esp_plots
subplot(1,2,2)
goods = find(~isnan(myData));
myData = griddata(lon(goods),lat(goods),myData(goods),lon,lat);
myData = myData.*mask;
clear mAx2
if strcmp(mp.plottype,'log')
    [mAx,mAx2,cbar] = logpcolorpsd(lon, lat, myData,mp.mrange,0);
elseif strcmp(mp.plottype,'linear')
    [mAx,cbar] = pcolorpsd(lon,lat,myData,mp.mrange,mp.tpoints);
end
cBarLabel = mp.title;
esp_plots
 print('-dpng',mp.print)
%  close