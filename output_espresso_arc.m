myFile = fullfile...
    ('C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\kN_pt5_zm5_bdir', ...
    'analysis1_gen_geo.mat');
myVars = {'mmed';'mIPR';'m95'};
fileNames = {'MAB_median';'MAB_hIPR';'MAB_95th_perc'};


myData = load(myFile,'lon','lat');
lon = myData.lon;
lat = myData.lat;

%Set up the grid
x = lon;
y = lat;
[m n] = size(x);
x = 0.5*(x(:,1:n-1) + x(:,2:n));
y = 0.5*(y(:,1:n-1) + y(:,2:n));
x = 0.5*(x(1:m-1,:) + x(2:m,:));
y = 0.5*(y(1:m-1,:) + y(2:m,:));

x1 = x(:,1) - (x(:,2)-x(:,1));
xE = x(:,end) + (x(:,end)-x(:,end-1));
x = [x1 x xE];
x1 = x(1,:) - (x(2,:) - x(1,:));
xE = x(end,:) + (x(end,:) - x(end-1,:));
x = [x1; x; xE];

y1 = y(:,1) - (y(:,2)-y(:,1));
yE = y(:,end) + (y(:,end)-y(:,end-1));
y = [y1 y yE];
y1 = y(1,:) - (y(2,:) - y(1,:));
yE = y(end,:) + (y(end,:) - y(end-1,:));
y = [y1; y; yE];

cd('C:\Users\sdalyander\Documents\StressAnalysis\Arc_MAB')
for vv = 1:length(myVars)
    thisVar = myVars{vv};
    if strcmp(thisVar,'mIPR')
        load(myFile,'m84','m16')
        myData = 0.5*(m84-m16); clear m84 m16
    else
    myData = load(myFile,thisVar);
    eval(['myData = myData.' thisVar ';'])
    end
    myData = squeeze(myData(1,:,:,:));
    
    dataSlice = squeeze(myData(1,:,:));
    goods = ~isnan(dataSlice);
    
    patches_x = NaN(5,numel(lon));
    patches_y = NaN(5,numel(lon));

    for ww = 1:numel(dataSlice)
        if isnan(dataSlice(ww)), continue,end
        [ii_w,jj_w] = ind2sub(size(dataSlice),ww);
        patches_x(:,ww) = [x(ii_w,jj_w) x(ii_w+1,jj_w) x(ii_w+1,jj_w+1) x(ii_w,jj_w+1) x(ii_w,jj_w)];
        patches_y(:,ww) = [y(ii_w,jj_w) y(ii_w+1,jj_w) y(ii_w+1,jj_w+1) y(ii_w,jj_w+1) y(ii_w,jj_w)];
    end
    
    clear Stress dataSlice
    
    dateNew = NaN(numel(goods),5);
    
    for mm = 1:5
        dataSlice = squeeze(myData(mm,:,:));
        dataSlice = dataSlice(:);
        dataSlice = dataSlice(goods);
        dataSlice(isnan(dataSlice)) = -9999; %For Esri
        dataNew(:,mm) = dataSlice; clear dataSlice
    end
    patches_x = patches_x(:,goods);
    patches_y = patches_y(:,goods);
    clear goods mm ii_w jj_w 
    
    for ww = 1:size(dataNew,1)
        Stress(ww).Geometry = 'Polygon';
        Stress(ww).Year = dataNew(ww,1);
        Stress(ww).Winter = dataNew(ww,2);
        Stress(ww).Spring = dataNew(ww,3);
        Stress(ww).Summer = dataNew(ww,4);
        Stress(ww).Fall = dataNew(ww,5);
        Stress(ww).Lon = patches_x(:,ww);
        Stress(ww).Lat = patches_y(:,ww);
    end
    shapewrite(Stress,[fileNames{vv} '.shp'])
end
    
%% Points
ncclear
load(fullfile('C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\point_stress', ...
    'compare_sediment_percentiles_to_skinstress_tseries_cohesive'), 'locs','obs_model_percentiles')
clear Stress
cc = 1;
obs_model_percentiles(isnan(obs_model_percentiles)) = -9999;
for ww = 1:size(locs,1)
    if (obs_model_percentiles(ww,1) == -9999), continue, end
    Stress(cc).Geometry = 'Point';
    Stress(cc).Lon = double(locs(ww,1));
    Stress(cc).Lat = double(locs(ww,2));
    Stress(cc).Year = double(obs_model_percentiles(ww,1));
    Stress(cc).Winter = double(obs_model_percentiles(ww,2));
    Stress(cc).Spring = double(obs_model_percentiles(ww,3));
    Stress(cc).Summer = double(obs_model_percentiles(ww,4));
    Stress(cc).Fall = double(obs_model_percentiles(ww,5));
    cc = cc+1;
end
shapewrite(Stress,'MAB_mobile_perc.shp')

%%
ncclear
load(fullfile('C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\point_stress', ...
    'mobility_event_analysis'), 'locs','mean_rec_interval')
clear Stress
cc = 1;
mean_rec_interval(isinf(mean_rec_interval)) = -9999;
mean_rec_interval(isnan(mean_rec_interval)) = -9999;
for ww = 1:size(locs,1)
    if (mean_rec_interval(ww,1) == -9999), continue, end
    Stress(cc).Geometry = 'Point';
    Stress(cc).Lon = double(locs(ww,1));
    Stress(cc).Lat = double(locs(ww,2));
    Stress(cc).Year = double(mean_rec_interval(ww,1));
    Stress(cc).Winter = double(mean_rec_interval(ww,2));
    Stress(cc).Spring = double(mean_rec_interval(ww,3));
    Stress(cc).Summer = double(mean_rec_interval(ww,4));
    Stress(cc).Fall = double(mean_rec_interval(ww,5));
    cc = cc+1;
end
shapewrite(Stress,'MAB_mobile_freq.shp')

    
