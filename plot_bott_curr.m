%%plot_mean_curr

ncclear, clc
load bott_current

%% Analysis
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';

tic
nc = cfdataset(mp.urlH);
var = nc.variable('tauwc');

count = zeros(size(lon));
msum = zeros(size(lon));


thresh = 0.5;
for tt = 1:length(time)
    disp(['On tt = ' num2str(tt) ' of ' num2str(length(time))])
    tauwc = squeeze(double(var.data(tt,:,:)));
    test = tauwc >= thresh;
    test = double(test);
    count = count + test;
    msum = msum+squeeze(mean_curr(tt,:,:)).*test;
end
clear var nc
toc

zVar = msum./count;


%% Load
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

%Pull grain size data
gFile = fullfile('C:\Users\sdalyander\Documents\SedimentTexture\FromBrad', ...
    'ecstdb2005.xls');
[nums,tex] = xlsread(gFile);
mc = nums(:,strcmp(tex(1,:),'MEAN'));
my = nums(:,strcmp(tex(1,:),'LATITUDE'));
mx = nums(:,strcmp(tex(1,:),'LONGITUDE'));
load(fullfile('C:\Users\sdalyander\Documents\Presentations\Geohab2011\figures', ...
    'MyColormaps'),'sed_texture2');
col2 = sed_texture2;

%% plot
LI =  [-74.4450  -71.4986   39.5121   41.7748];
isIn = find(lon >= LI(1) & lon <= LI(2) & lat >= LI(3) & lat <= LI(4));
[ii2,jj2] = ind2sub(size(lon),isIn);
close all, figure
Umag = abs(zVar);
phideg = angle(zVar)*180/pi;
[mAx,cbar] = pcolorpsd(lon,lat,Umag,[0 0.1], [0:0.025:0.1]), hold on
% cBarLabel = 'Bottom Current (m/s), \tau_w_c > 0.5 Pa';
cBarLabel = 'Mean Bottom Current (m/s)';
esp_plots
set(mAx,'xlim',LI(1:2),'ylim',LI(3:4));
set(gca,'xlim',LI(1:2),'ylim',LI(3:4))
sp = 1;
sc = 1.5;
ii = min(ii2):max(ii2); jj = min(jj2):max(jj2);
lonCut = lon(ii,jj); latCut = lat(ii,jj);
UmagCut = Umag(ii,jj); phidegCut = phideg(ii,jj);
a = arrowsafe(lonCut,latCut,...
    sc*UmagCut,phidegCut,1.5e-2);
% a = arrowsafe(lon(1:sp:end,1:sp:end),lat(1:sp:end,1:sp:end),...
%     sc*Umag(1:sp:end,1:sp:end),phideg(1:sp:end,1:sp:end),1e-5);

arrowsafe
set(a,'color',[0.5 0.5 0.5])
% ax = -70.2; ay = 36.5;
ax = -74.20; ay = 41.54;
a1 = arrowsafe(ax,ay,sc*0.2,0,1.5e-2);
set(a1,'color','k')
t = text(ax+sc*0.2/2,ay-0.1,'20 cm/s');
set(t,'fontname','times','vert','top','horiz','center','fontsize',14)


%% Grain Size
hold on
cRange = [0 6];
mc2 = 1+round((mc - cRange(1))/(cRange(2) - cRange(1))*63);
mc2(mc2 < 1) = 1;
mc2(mc2 > 64) = 64;
bads = isnan(mc2);
mx(bads) = [];
my(bads) = [];
mc2(bads) = [];
mc(bads) = [];
big = find(mc <= 1);
med = find(mc >1 & mc < 4);
small = find(mc >= 4);
ms = 3;
for ii = 1:length(med)
    pp = med(ii);
    p(ii) = plot(mx(pp),my(pp),'o');
    set(p(ii),'color',col2(mc2(pp),:),'markerface',col2(mc2(pp),:),'markersize',ms)
end
for ii = 1:length(big)
    pp = big(ii);
    p2(ii) = plot(mx(pp),my(pp),'o');
    set(p2(ii),'color',col2(mc2(pp),:),'markerface',col2(mc2(pp),:),'markersize',ms)
end
for ii = 1:length(small)
    pp = small(ii);
    p3(ii) = plot(mx(pp),my(pp),'o');
    set(p3(ii),'color',col2(mc2(pp),:),'markerface',col2(mc2(pp),:),'markersize',ms)
end
