%% mean_curr_stress_crit

ncclear
clc
mp.urlC = ...
   'http://geoport.whoi.edu/thredds/dodsC/coawst_2_2/fmrc/coawst_2_2_best.ncd';
% mp.urlC = ...
%    'http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2009_da/his';
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';

nc=cfdataset(mp.urlC);
% timeC = nc.time('ocean_time');
nc2 = mDataset(mp.urlC);
timeC = nj_time(nc2,'ubar');
close(nc2), clear nc2
lon = nc.data('lon_rho');
lat = nc.data('lat_rho');
h = nc.data('h');
mask = nc.data('mask_rho');

cin = find(timeC >= datenum(2010,05,01,00,00,00) & ...
    timeC <= datenum(2011,05,01,00,00,00));
timeC = timeC(cin);   %Cushion
LI =  [-74.4450  -71.4986   39.5121   41.7748];
ind =  find(lon(:) >= LI(1) & lon(:) <= LI(2) & lat(:) >= LI(3) & ...
    lat(:) <= LI(4));
[ii,jj] = ind2sub(size(lon),ind);
ii = min(ii):max(ii);
jj = min(jj):max(jj);

varU = nc.variable('u');
u = double(squeeze(varU.data(cin,1,ii,[min(jj)-1 jj])));
% varU = nc.variable('ubar');
% u = double(squeeze(varU.data(cin,ii,[min(jj)-1 jj])));
Ucur = u2rho_3d(u);
Ucur(:,:,1) = []; Ucur(:,:,end) = [];
% varV = nc.variable('vbar');
varV = nc.variable('v');
try
v = double(squeeze(varV.data(cin,1,[min(ii)-1 ii],jj)));
% v = double(squeeze(varV.data(cin,[min(ii)-1 ii],jj)));
Vcur = v2rho_3d(v);
Vcur(:,1,:) = []; Vcur(:,end,:) = [];
catch
%         v = double(squeeze(varV.data(cin,[min(ii)-1:1:81],jj)));
    v = double(squeeze(varV.data(cin,1,[min(ii)-1:1:81],jj)));
Vcur = v2rho_3d(v);
Vcur(:,1,:) = [];
end

var = nc.variable('angle');
angle_r = var.data(:);
clear u v var nc varU varV

mask = mask(ii,jj);
h = h(ii,jj);
lon = lon(ii,jj);
lat = lat(ii,jj);
angle_r = angle_r(ii,jj);

mean_u = squeeze(nanmean(Ucur,1));
mean_v = squeeze(nanmean(Vcur,1));
Umag = sqrt(mean_v.^2 + mean_u.^2);
phic = atan2(mean_v,mean_u) + angle_r;
clear mean_u mean_v
mean_curr = (Umag.*cos(phic)) + 1i.*(Umag.*sin(phic));
% clear Umag phic varV varU u v


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

mask2 = mask;
mask2(~mask2) = NaN;

%% plot
LI =  [-74.4450  -71.4986   39.5121   41.7748];
isIn = find(lon >= LI(1) & lon <= LI(2) & lat >= LI(3) & lat <= LI(4));
[ii2,jj2] = ind2sub(size(lon),isIn);
close all, figure
Umag = abs(mean_curr).*mask2;
phideg = angle(mean_curr)*180/pi;
[mAx,cbar] = pcolorpsd(lon,lat,Umag,[0 0.1], [0:0.025:0.1]), hold on
cBarLabel = 'Mean Bottom Current, COAWST (m/s)';
% cBarLabel = 'Mean Depth-Averaged Current, ESPreSSO (m/s)';
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



