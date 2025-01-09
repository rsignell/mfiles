ncclear
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';

nc = mDataset(mp.urlH);
time = nj_time(nc,'tauwc');
lon = nc{'lon'}(:);
lat = nc{'lat'}(:);

test = squeeze(double(nc{'tauwc'}(100,:,:)));

p = [-69.8837   41.0680];
ind = nearxy(lon(:),lat(:),p(1),p(2));
[jj,ii] = ind2sub(size(lon),ind);

tauc = nc{'tauc'}(:,jj,ii);

[tauc_lp,time_lp] = plfilt(tauc,time,6);

%% Plot
close all
figure
subplot(2,1,1)
plot(time,tauc,'k')
hold on
plot(time_lp,tauc_lp,'r')
xlim([datenum(2011,01,01) datenum(2011,02,01)])