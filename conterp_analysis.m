%conterp_analysis

ncclear
%% Grain size
gFile = fullfile('C:\Users\sdalyander\Documents\SedimentTexture\FromBrad', ...
    'ecstdb2005.xls');
[nums,tex] = xlsread(gFile);
meansed = nums(:,strcmp(tex(1,:),'MEAN'));
my = nums(:,strcmp(tex(1,:),'LATITUDE'));
mx = nums(:,strcmp(tex(1,:),'LONGITUDE'));
p = [mx(:) my(:)]; clear mx my nums tex
 load C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\perc_grainsize_crit_smaller_percentiles box_obs_in
 isIn = ones(size(p,1),1);
isIn(isnan(box_obs_in(:,1))) = NaN;

%Pick the file, variable to look at
mp.file = 'analysis1_gen_geo.mat';
mp.pvar = 'mmean(1,1,:,:)';   %(ttype,timeperiod,lon,lat)
% mp.pvar = 'm95(1,1,:,:)';
mp.urlC = ...
   'http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2009_da/his';
mp.depthVar = 'h';

nc = cfdataset(mp.urlC);
var = nc.variable(mp.depthVar);
bathy = var.data(:);
var = nc.variable('lon_rho');
lon = var.data(:);
var = nc.variable('lat_rho');
lat = var.data(:);
var = nc.variable('mask_rho');
mask = var.data(:);
clear nc var

mask(~mask) = NaN;
mask(bathy > 375) = NaN;
load(mp.file);
eval(['myData = squeeze(' mp.pvar ');'])

%% do it
close all
figure(100), hold on
subplot(3,1,1)
pos1 = get(gca,'position');
subplot(3,1,2)
pos2 = get(gca,'position');
subplot(3,1,3)
pos3 = get(gca,'position');
cList = {'b'; 'c'; 'g'; 'r'; 'm'};
ccList = [60:10:100];
clear p2
goods = find(~isnan(meansed));
fSed = TriScatteredInterp(p(goods,1),p(goods,2),meansed(goods));
for cci = 1:length(ccList)
cc = ccList(cci);
figure(1)
[xi,yi,wi] = conterp_psd(lon,lat,bathy,cc,myData);
figure(100)
subplot(3,1,1), hold on
bads = find(xi > -60 | yi > 41.5);
wi(bads) = NaN;
p2(cci) = plot(xi,wi,cList{cci});
labs{cci} = [num2str(ccList(cci)) ' m'];
subplot(3,1,2), hold on
plot(xi(1:end-1),diff(wi),cList{cci});
figure(100)
subplot(3,1,3), hold on
si = fSed(xi,yi);
si(bads) = NaN;
p3(cci) = plot(xi,si,cList{cci});
isGood = find(~isnan(si) & ~isnan(wi));
[r4,p4] = corrcoef(wi(isGood),si(isGood))
end
figure(100), subplot(3,1,1)
l = legend(p2,labs,'location','northwest');
set(l,'fontsize',6,'box','off','color','none')
xlim([-74 -68.5])
% view(180,-90)
xlabel('<----West      East------>')
title('Mean')
box on
ylim([0 0.4])
ylabel('\tau_w_c (Pa)')
figure(100), subplot(3,1,2)
% l = legend(p2,labs,'location','northwest');
% set(l,'fontsize',6,'box','off','color','none')
xlim([-74 -68.5])
% view(180,-90)
xlabel('<----West      East------>')
title('\DeltaMean')
box on
ylim([-0.01 0.01])
ylabel('\Delta\tau_w_c (Pa)')
    figure(100),subplot(3,1,3)
% l = legend(p3,labs,'location','eastoutside');
xlim([-74 -68.5])
% view(180,-90)
box on
ylabel('Sediment Texture (\phi)')
l = line([min(xlim) max(xlim)], [4 4]);
set(l,'color','r')
xlabel('<----West      East------>')
title('Sediment Texture (Mean)')
ylim([-1 8])
close(1)


figure, hold on
pcolorjw(lon,lat,bathy.*mask)
[c1,h1] = contour(lon,lat,bathy,[10:10:380]);
set(h1,'color','w')
for cci = 1:length(ccList)
[c,h] = contour(lon,lat,bathy,[ccList(cci) ccList(cci)]);
set(h,'color',cList{cci},'linewidth',1)
end
caxis([0 375])
colormap('gray')
colorbar
xlim([-74 -68.5])
ylim([38 42])
dasp(41)
legend(h1,'10 m contours')
