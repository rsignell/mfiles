%% Plot results
clear
mp.dir = 'C:\Users\sdalyander\Documents\StressAnalysis\SABGOM\stress_analysis';
mp.pdir = 'C:\Users\sdalyander\Documents\StressAnalysis\SABGOM\stress_analysis';
cd(mp.dir)

%% Mean
mp.file = 'analysis1_gen_geo.mat';
mp.pvar = '(m84(5,4,:,:)-m16(5,4,:,:))./2'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
mp.title = 'hIPR Non-Tidal Current Stress, Summer'; %<--------CHANGE THIS
mp.plottype = 'log';    %Log or linear
mp.plotclean = 'sab_plots';
mp.mrange = [0.01 2];
doCut = 0;
mp.tpoints = [0:0.5:2];
mp.print = 'mab_sum_taut_med'; %<--------CHANGE THIS

%% Load up everything
%Data
cd(mp.dir)
load(fullfile(mp.dir,mp.file))
load(fullfile(mp.dir,'mask'))

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


nc = mDataset('dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CFSR_7grid_v2/swan_agg.ncml');
lonAll = nc{'lon'}(:);
latAll = nc{'lat'}(:);
depthAll = nc{'depth'}(:);
close(nc), clear nc

%% Plots
%Feed in data
close all
eval(['myData = mask.*squeeze(' mp.pvar ');'])
% eval(['myData = squeeze(' mp.pvar ');'])
myData(myData == Inf) = NaN;    %Recurrence interval
% For Gulf of Mexico plots
%  [jj,ii] = lonlat2ij(lon,lat,[-98 -80 22 32]);
%  Lon = lon(jj,ii);
%  Lat = lat(jj,ii);
%  myData = myData(jj,ii);
% For South Atlantic Bight plots
 [jj,ii] = lonlat2ij(lon,lat,[-85 -70 20 38]);
 Lon = lon(jj,ii);
 Lat = lat(jj,ii);
 myData = myData(jj,ii);
figure, hold on
clear mAx2
if strcmp(mp.plottype,'log')
    [mAx,mAx2,cbar] = logpcolorpsd(Lon, Lat, myData,mp.mrange,0);
elseif strcmp(mp.plottype,'linear')
    [mAx,cbar] = pcolorpsd(Lon,Lat,myData,mp.mrange,mp.tpoints,doCut);
end
cBarLabel = mp.title;
eval(mp.plotclean)
box off