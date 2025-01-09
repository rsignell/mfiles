%% Plot results
ncclear
mp.dir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\point_stress';
mp.pdir = 'C:\Users\sdalyander\Documents\MyPublications\MAB_Paper\Figures';
cd(mp.dir)

load ../kN_pt5_zm5_bdir/analysis1_gen_geo lon lat
% 
% %% Magnitude of stress
% mp.file = 'polygons\polyout';
% mp.pvar = 'magStress'; %<------CHANGE THIS 
% mp.title = 'Stress Magnitude (95th Percentile)'; %<--------CHANGE THIS
% mp.plotclean = 'esp_plots';
% mp.print = 'magStress'; %<--------CHANGE THIS
% mp.list = {'Region:'; 'Low (< 0.5 Pa)'; 'Intermediate (0.5-1 Pa)'; 'High (> 1 Pa)'};
% doCut = 1;
% 
%% Bed mobility
mp.file = 'polyout';
mp.pvar = 'freqStress'; %<------CHANGE THIS 
mp.title = 'Stress Frequency (Ratio of Low/High Frequency Events)'; %<--------CHANGE THIS
mp.plotclean = 'esp_plots';
mp.print = 'freqStress'; %<--------CHANGE THIS
mp.list = {'Region:'; 'Infrequent (> 10)'; 'Intermediate (2-10)'; 'Frequent (< 2)'};

% % Bed mobility
% mp.file = 'polygons\polyout';
% mp.pvar = 'magDisturb'; %<------CHANGE THIS 
% mp.title = 'Mobility Magnitude (90th Percentile Stress)'; %<--------CHANGE THIS
% mp.plotclean = 'esp_plots';
% mp.print = 'magDisturb'; %<--------CHANGE THIS
% mp.list = {'Region:'; 'Low (< 25%)'; 'Intermediate (25-75%)'; 'High (> 75%)'};
% % 
% %% Bed mobility
% mp.file = 'polyout';
% mp.pvar = 'freqDisturb'; %<------CHANGE THIS 
% mp.title = 'Mobility Frequency (20th Percentile Bed)'; %<--------CHANGE THIS
% mp.plotclean = 'esp_plots';
% mp.print = 'freqDisturb'; %<--------CHANGE THIS
% mp.list = {'Region:'; 'Infrequent (< 5%)'; 'Intermediate (5-50%)'; 'Frequent (> 50%)'};

%% Load up everything
%Data
cd(mp.dir)
load(fullfile(mp.dir,mp.file))
load ../kN_pt5_zm5_bdir/mask

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
mp.mrange = [1 3];
mp.tpoints = [1:3];
eval(['myData = squeeze(' mp.pvar ').*mask;'])
isOne = (myData == 1);
isThree = (myData == 3);
myData(isOne) = 3;
myData(isThree) = 1;
myData(myData == Inf) = NaN;    %Recurrence interval
close all
figure
clear mAx2
[mAx,cbar] = pcolorpsd(lon,lat,myData,mp.mrange,mp.tpoints,0);
cBarLabel = mp.title;


esp_plots_poly
%  print('-dpng',fullfile(mp.pdir,mp.print))
  print('-depsc2','-tiff',fullfile(mp.pdir,mp.print))
%  close