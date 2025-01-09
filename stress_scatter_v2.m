%stress_scatter_v2

ncclear, close all, clc

gp = 90;
sp = 10;

eval(['load grid_sed_texture_' num2str(gp) 'p_thresh'])
load analysis1_gen_geo

%% Load up everything
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
close all
if sp ~= 50
eval(['xVar = squeeze(m' num2str(sp) '(1,1,:,:));'])
else
    xVar = squeeze(mmed(1,1,:,:));
end
eval(['xlab = ''' num2str(sp) 'th Percentile Stress (Pa)'';']);
eval(['yVar = crit_' num2str(gp) 'p_stress;'])
eval(['ylab = ''\tau_c_r_i_t for ' num2str(gp) 'th Percentile Grain Size (Pa)'';']);

figure, hold on
bads = find(yVar < xVar);
plot(xVar(:),yVar(:),'k.')
plot(xVar(bads),yVar(bads),'mo')
axis equal
axis([0 2 0 2])
l = line([0 2], [0 2]);
set(l,'color','r')
xlabel(xlab)
ylabel(ylab)
box on


figure, hold on
mp.mrange = [1e-2 4];
clear mAx2
[mAx,mAx2,cbar] = logpcolorpsd(lon, lat,xVar,mp.mrange,0);
cBarLabel = xlab;
esp_plots
p = plot(lon(bads),lat(bads),'o');
set(p,'color','w','markerfacecolor','m','markersize',5)

%% print
figure(1), eval(['print -dpng scat_crit' num2str(gp) 'p_' num2str(sp) 'thstress;']); % close
figure(2), eval(['print -dpng s' num2str(sp) 'th_outlier_crit' num2str(gp) 'p;']);  %close

