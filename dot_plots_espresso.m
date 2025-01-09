%grain_size_plot
ncclear
mp.pdir = 'C:\Users\sdalyander\Documents\MyPublications\MAB_Paper\Figures';
load ../kN_pt5_zm5_bdir/analysis1_gen_geo.mat lon lat 
% % 
%% fix the stress percentile, look at the sediment
%  load perc_grainsize_crit_smaller_percentiles
% thresh = 80;
% cRange = [0 100];
% mc = perc_crit_smaller_obs_all(:,thresh/10);
% mTitle = ['Fraction of Bed Mobile ' num2str(thresh) '% of the Time'];
% mTicks = 0:25:100;
% 
% fix the sediment percentile, look at the stress
load compare_sediment_percentiles_to_skinstress_tseries_cohesive
cRange = [0 2];
mc = log10(obs_model_percentiles(:,4));
mTitle = ['Summer Percentage of Time Bed is Mobilized'];
mTicks = [-1:1:2];
mp.print = 'perc_mobile_sum.eps';
doFlip = 0;
p = locs;

% %% fix the sediment percentile, look at the stress
% load mobility_frequency_cohesive
% % thresh = 80;
% cRange = [0 100];
% mc = 100*tide_frac_all;
% mTitle = ['Tide Frac, All'];
% % mTitle = ['Tide Frac for >' num2str(thresh) '% of Bed is Mobilized'];
% mTicks = 0:25:100;

% % % Grain size distribution
% % Pull grain size data
% % gFile = fullfile('C:\Users\sdalyander\Documents\SedimentTexture\FromBrad', ...
% %     'ecstdb2005.xls');
% % [nums,tex] = xlsread(gFile);
% % stdsed = nums(:,strcmp(tex(1,:),'STDEV'));
% % meansed = nums(:,strcmp(tex(1,:),'MEAN'));
% % pct_grav = nums(:,strcmp(tex(1,:),'GRAVEL_PCT'));
% % pct_sand = nums(:,strcmp(tex(1,:),'SAND_PCT'));
% % pct_silt = nums(:,strcmp(tex(1,:),'SILT_PCT'));
% % pct_clay = nums(:,strcmp(tex(1,:),'CLAY_PCT'));
% % stdsed = (2.^(-1*stdsed))/1000;
% % my = nums(:,strcmp(tex(1,:),'LATITUDE'));
% % mx = nums(:,strcmp(tex(1,:),'LONGITUDE'));
% % p = [mx(:) my(:)]; clear mx my nums tex
% % isIn = ones(size(p,1),1);
% % isIn(isnan(box_obs_in(:,1))) = NaN;

%%
% cRange = [1 5];
% mTicks = [1:1:5];
% mc = meansed.*isIn;
% mTitle = 'Mean Grain Size (\phi)';

%% process contributions
% load mobility_process
% thresh = 20;
% cRange = [0 1];
% mc = res_mean(:,thresh/10)./total_mean(:,thresh/10);
% mTitle = ['Wave Fraction of Mobility'];
% % mTitle = ['Tide Frac for >' num2str(thresh) '% of Bed is Mobilized'];
% mTicks = 0:0.25:1;

% %% mobility frequency
% load mobility_frequency
% thresh = 20;
% cRange = [0 30];
% mc = squeeze(spec_analysis_low2high);
% mTitle = ['Low/High'];
% % mTitle = ['Tide Frac for >' num2str(thresh) '% of Bed is Mobilized'];
% mTicks = 0:10:30;

% %% mobility event analysis
% load mobility_event_analysis
% % thresh = 20;
% month = 4;
% cRange = [0 10];
% mc = squeeze(mean_rec_interval(:,month));
% mc(mc == Inf) = NaN;
% mTitle = ['Summer Mean Recurrence Interval (Days)'];
% mTicks = 0:2:10;
% mp.print = 'summ_recurr_int.eps';
% doFlip = 1;
% p = locs;

% %% mobility event analysis
% load iscohesive
% % thresh = 20;
% % month = 4;
% cRange = [0 1];
% mc = isCo;
% mc(isnan(Cboxes(:,1)))= NaN;
% mc(isnan(csum(:,end))) = NaN;
% % mc(mc == Inf) = NaN;
% mTitle = 'Cohesive';
% mTicks = [0 1];
% mp.print = 'isCo.eps';
% doFlip = 0;
% p = double(locs);
% tau_crit = log(tau_crit);

% mobility event analysis
% load mobility_event_analysis_daily
% thresh = 80;
% month = 4;
% cRange = [0 10];
% mc = squeeze(avg_days_between_events(:,month,thresh/10));
% mc(mc == Inf) = NaN;
% mTitle = ['Summer Days Betweem Events'];
% mTicks = 0:2:10;



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

%States
PBFile = fullfile('C:\Users\sdalyander\Documents\Maps\WDB','namer-pby.txt');
PB = read_wdbfile(PBFile);
PBFile = fullfile('C:\Users\sdalyander\Documents\Maps\WDB','namer-cil.txt');
PB2 = read_wdbfile(PBFile);

% cd(mp.mdir)

xfac=cos(38*pi/180);

%% Visualize
% close all
figure, hold on
col2 = colormap;
if doFlip
col2 = flipud(col2);
% col2(end,:) = 0;
else
%     col2(end,:) = 0;
end
set(gcf,'colormap',col2)
ax = [min(min(lon)) max(max(lon)) min(min(lat)) max(max(lat))];
axis(ax)
mAx = gca;
caxis(cRange)
cbar = colorbar('peer',mAx,'location','southoutside');
set(cbar,'xtick',mTicks);
% caxis([0 6])
cBarLabel = mTitle;
esp_plots

hold on
mc2 = 1+round((mc - cRange(1))/(cRange(2) - cRange(1))*63);
mc2(mc2 < 1) = 1;
mc2(mc2 > 64) = 64;
% bads = isnan(mc2);
% p(bads,:);
% mc2(bads) = [];
% mc(bads) = [];

for ii = 1:length(mc2)
    if isnan(mc2(ii)), continue, end
    p2(ii) = plot(p(ii,1),p(ii,2),'o');
    set(p2(ii),'color',col2(mc2(ii),:),'markerface',col2(mc2(ii),:),'markersize',4)
end
  print('-depsc2','-tiff',fullfile(mp.pdir,mp.print))
%  close