%% Plot results
ncclear
mp.dir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\kN_pt5_zm5_bdir';
% mp.dir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\kN_pt5_zm5';
mp.pdir = 'C:\Users\sdalyander\Documents\MyPublications\MAB_Paper\Figures';
cd(mp.dir)

% %% Mean
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = 'mmed(4,4,:,:)'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
% mp.title = 'Summer Median Tidal Stress (Pa)'; %<--------CHANGE THIS
% mp.plottype = 'log';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [1e-4 1];
% mp.print = 'mab_sum_taut_med'; %<--------CHANGE THIS

% %% Wave/Current fractions
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = 'taur_percc(end,1,:,:)'; %<------CHANGE THIS (gsize,timeperiod,lon,lat) 
% % mp.title = 'Yearly Mean \tau_c_,_T_i_d_e/\tau_c'; %<--------CHANGE THIS
% mp.title = 'Yearly Mean \tau_c_,_R_e_s/\tau_c'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 10];
% mp.tpoints = [0:2.5:10];
% mp.print = 'mab_yrly_tfrac_62micron.png'; %<--------CHANGE THIS

% % Wave/Current fractions
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = 'taut2s(end,1,:,:)'; %<------CHANGE THIS (gsize,timeperiod,lon,lat) 
% mp.title = 'Yearly Mean \tau_c_,_T_i_d_e/\tau_c'; %<--------CHANGE THIS
% mp.title = 'Yearly Mean \tau_c_,_T_i_d_e/(\tau_c_,_R_e_s+\tau_w)'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 10];
% mp.tpoints = [0:2.5:10];
% mp.print = 'mab_yrly_tfrac_62micron.png'; %<--------CHANGE THIS
% 
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = 'tauw_perc(end,1,:,:)'; %<------CHANGE THIS (gsize,timeperiod,lon,lat) 
% mp.title = 'Yearly Mean \tau_w_/\tau_w_c'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 1];
% mp.tpoints = [0:0.25:1];
% mp.print = 'mab_yrly_tfrac_62micron.png'; %<--------CHANGE THIS

% %% STD
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = '(m84(1,1,:,:)-m16(1,1,:,:))/2'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
% mp.title = 'Yearly St. Deviation \tau_w_c (Pa)'; %<--------CHANGE THIS
% mp.plottype = 'log';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [1e-2 4];
% mp.print = 'mab_yrly_spread'; %<--------CHANGE THIS


% %% STD
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = 'mstd(1,1,:,:)'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
% mp.title = 'Yearly St. Deviation \tau_w_c (Pa)'; %<--------CHANGE THIS
% mp.plottype = 'log';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [1e-2 4];
% mp.print = 'mab_yrly_std'; %<--------CHANGE THIS

%% 95%
mp.file = 'analysis1_gen_geo.mat';
mp.pvar = 'm95(4,4,:,:)'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
mp.title = 'Summer 95^t^h Percentile \tau_t (Pa)'; %<--------CHANGE THIS
mp.plottype = 'log';    %Log or linear
mp.plotclean = 'esp_plots';
mp.mrange = [1e-1 100];
mp.print = 'mab_sum_p95_taut.eps'; %<--------CHANGE THIS
% 
% %% Cross correlation
% mp.file = 'analysis3_corr_stress2NAM';
% mp.pvar = 'xcor_wind_wave(1,:,:)';  %<-----CHANGE THIS (mont,lon,lat)
% mp.title = 'Correlation of Wave Stress to Local Wind Stress, All Year'; %<-----CHANGE THIS
% mp.plottype = 'linear';
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 0.8];
% mp.tpoints = [0:0.2:0.8];
% doCut = 0;
% mp.print = 'mab_xcor_wind_wave_yrly';    %<-------CHANGE THIS

% %% Spectral analysis
% mp.file = 'analysis3_spec_welch.mat';
% mp.pvar = 'spec_analysis_low(1,1,:,:)./spec_analysis_high(1,1,:,:)';  %<-----CHANGE THIS (ttype,month,lon,lat)
% mp.title = 'Ratio of Low to High Frequency Energy for twc'; %<-----CHANGE THIS
% mp.plottype = 'log';
% mp.plotclean = 'esp_plots';
% mp.mrange = [10^-1 75];
% mp.tpoints = [0:2:10];
% mp.print = 'mab_spec_tide_low2high';    %<-------CHANGE THIS

%% Bed mobility
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = 'mbed_mobil(2,1,:,:)'; %<------CHANGE THIS (gsize,timeperiod,lon,lat) 
% mp.title = 'Yearly Bed Mobility, 1 mm'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 5];
% mp.tpoints = [0:1:5];
% % mp.print = 'mab_yrly_bm_1mm'; %<--------CHANGE THIS
% mp.print = 'test';

%% Stress Fractions
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = 'm95(2,2,:,:)./(m95(2,2,:,:)+m95(4,2,:,:)+m95(5,2,:,:))'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
% mp.title = 'Winter Mean \tau_w/Total Mean \tau_w_c'; %<--------CHANGE THIS
% % mp.title = 'Yearly Mean \tau_w/ Yearly Mean \tau_w_c'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 1];
% mp.tpoints = [0:0.25:1];
% doCut = 1;
% mp.print = 'mab_wint_wave_over_total'; %<--------CHANGE THIS
% 
% %% CV
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = '(m84(1,1,:,:) - m16(1,1,:,:))./(2*mmed(1,1,:,:))'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
% mp.title = 'Coefficient of Variation (unitless)'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 2];
% mp.tpoints = [0:0.5:2];
% mp.print = 'mab_yrly_tauwc_cv'; %<--------CHANGE THIS
% doCut = 0;

% %% GCV
% mp.file = 'analysis1_gen_geo.mat';
% mp.pvar = 'exp(log(gstd(1,1,:,:))./gmean(1,1,:,:))'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
% mp.title = 'Coefficient of Variation (unitless)'; %<--------CHANGE THIS
% mp.plottype = 'linear';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [0 2];
% mp.tpoints = [0:0.5:2];
% mp.print = 'mab_yrly_tauwc_cv'; %<--------CHANGE THIS
% doCut = 0;

% %% CV
% mp.file = 'analysis7_lowpass_highfreq';
% mp.pvar = 'lp_std(1,1,:,:)./hf_std(1,1,:,:)'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
% mp.title = 'Low/High Frequency'; %<--------CHANGE THIS
% mp.plottype = 'log';    %Log or linear
% mp.plotclean = 'esp_plots';
% % mp.mrange = [0 2];
% mp.tpoints = [10^-2 10^2];
% mp.print = 'mab_yrly_tauwc_hf_lp'; %<--------CHANGE THIS

%% Ellipticity
% mp.file = 'ellipticity';
% mp.pvar = 'abs(bot_elliptic)'; %<------CHANGE THIS (ttype,timeperiod,lon,lat) 
% mp.title = 'Ellipticity of M2 Tide'; %<--------CHANGE THIS
% mp.plottype = 'log';    %Log or linear
% mp.plotclean = 'esp_plots';
% mp.mrange = [10^0 10^2];
% mp.tpoints = [10^0 10^2];
% mp.print = 'elliptic'; %<--------CHANGE THIS

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



%% Plots
%Feed in data
eval(['myData = mask.*squeeze(' mp.pvar ');'])
myData(myData == Inf) = NaN;    %Recurrence interval
% close all
figure
clear mAx2
if strcmp(mp.plottype,'log')
    [mAx,mAx2,cbar] = logpcolorpsd(lon, lat, myData,mp.mrange,0);
elseif strcmp(mp.plottype,'linear')
    [mAx,cbar] = pcolorpsd(lon,lat,myData,mp.mrange,mp.tpoints,doCut);
end
cBarLabel = mp.title;
esp_plots
box off

    
%  print('-dpng',fullfile(mp.pdir,mp.print))
  print('-depsc2','-tiff',fullfile(mp.pdir,mp.print))
%  close