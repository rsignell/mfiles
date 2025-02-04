%%mean_currents

% tic
ncclear
clc
mp.urlC = ...
   'http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2009_da/his';
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';

nc=cfdataset(mp.urlC);
time = nc.time('ocean_time');
lon = nc.data('lon_rho');
lat = nc.data('lat_rho');

cin = find(time >= datenum(2010,05,01,00,00,00) & ...
    time <= datenum(2011,05,01,00,00,00));
time = time(cin);

%%Load


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

%Analysis:  mean depth-averaged current
var = nc.variable('ubar');
ubar = var.data(cin,:,:);
var = nc.variable('vbar');
vbar = var.data(cin,:,:);
var = nc.variable('angle');
angle = var.data(:);
clear nc, var

ubar = double(ubar);
vbar = double(vbar);

%% Analysis
mean_ubar = squeeze(nanmean(ubar,1));
mean_vbar = squeeze(nanmean(vbar,1));
mean_ubar = u2rho_2d(mean_ubar);
mean_vbar = v2rho_2d(mean_vbar);
Umag = sqrt(mean_ubar.^2 + mean_vbar.^2);
phic = atan2(mean_vbar,mean_ubar) + angle;
phideg = phic*180/pi;

nc = cfdataset(mp.urlH);
timeS = nc.time('time');


%% plot
close all, figure
[mAx,cbar] = pcolorpsd(lon,lat,Umag,[0 0.5], [0:0.1:0.5]), hold on
cBarLabel = 'Mean Current (m/s)';
esp_plots
sp = 4
sc = 4;
a = arrowsafe(lon(1:sp:end,1:sp:end),lat(1:sp:end,1:sp:end),...
    sc*Umag(1:sp:end,1:sp:end),phideg(1:sp:end,1:sp:end),1e-1);
%axis([min(lon(:)) max(lon(:)) min(lat(:)) max(lat(:))]), arrowsafe
set(a,'color',[0.7 0.7 0.7])
