%model_vals_at_sed_cores

%Run identify_cirt_stress_from_obs first

ncclear, clc
load analysis1_gen_geo
load perc_grainsize_crit_smaller_percentiles box_obs_in mc gsizel p

sList = [10:10:90];
gList = [10:10:90];

% Open the dataset
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';
nc = ncdataset(mp.urlH);
lon = nc.data('lon');
lat = nc.data('lat');
clear nc

model_percentile = NaN([size(mc,1) size(sList)]);
obs_percentiles = NaN([size(mc,1) size(gList)]);
obs_crit_perc = NaN([size(mc,1) size(gList)]);

for gg = 1:size(mc,1)
    ii = box_obs_in(gg,1);
    jj = box_obs_in(gg,2);
    
    if isnan(ii), continue, end

    for s = 1:length(sList)
        if sList(s) == 50
            model_percentile(gg,s) = squeeze(mmed(1,1,ii,jj));
        else
            eval(['model_percentile(gg,s) = m' num2str(sList(s),'%02.0f') ...
                '(1,1,ii,jj);'])
        end
    end
    
    gall = mc(gg,:); gall(isnan(gall)) = 0;
    if all(gall==0), continue, end
    gall = cumsum(gall);
        
    for g = 1:length(gList)
        tt = find(gall >= gList(g),1,'first');
        obs_percentiles(gg,g) = gsizel(tt);
        [~,obs_crit_perc(gg,g)] = pmsoulsby(obs_percentiles(gg,g),1);
    end
end

clear m05 m10 m20 m30 m40 m60 m70 m80 m90 sum_mc
clear mmed mbed_mobil mmean mstd mlists gg ii jj ind isIn isOn
clear bads first g gall last m m25 m75 m95 mc mspread n nodes s
clear tauc_perc taur_perc taut_perc tauw_perc time tt ttypes x y 

model_percentile = squeeze(model_percentile);
obs_crit_perc = squeeze(obs_crit_perc);
obs_percentiles = squeeze(obs_percentiles);

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
sp = 50;
gp = 90;
if sp == 50
    load analysis1_gen_geo mmed
    sVar = squeeze(mmed(1,1,:,:));
else
eval(['load analysis1_gen_geo m' num2str(sp,'%02.0f')])
eval(['sVar = squeeze(m' num2str(sp,'%02.0f') '(1,1,:,:));'])
end
xVar = model_percentile(:,sp/10);
yVar = obs_crit_perc(:,gp/10);
eval(['xlab = ''' num2str(sp) 'th Percentile Stress (Pa)'';']);
eval(['ylab = ''\tau_c_r_i_t for ' num2str(gp) 'th Percentile Grain Size (Pa)'';']);

figure, hold on
bads = find(yVar < xVar);
goods = find(yVar >= xVar);
plot(xVar(:),yVar(:),'k.')
plot(xVar(bads),yVar(bads),'mo')
axis equal
axis([0 4 0 4])
l = line([0 4], [0 4]);
set(l,'color','r')
xlabel(xlab)
ylabel(ylab)
box on

figure, hold on
mp.mrange = [1e-2 2];
clear mAx2
[mAx,mAx2,cbar] = logpcolorpsd(lon, lat,sVar,mp.mrange,0);
cBarLabel = ['Stress Exceeded ' num2str(100-sp,'%.0f') '% of the Time'];
esp_plots
p2 = plot(p(bads,1),p(bads,2),'o');
set(p2,'color','m','markerfacecolor','m','markersize',5)
p3 = plot(p(goods,1),p(goods,2),'o');
set(p3,'color','w','markersize',5)
        