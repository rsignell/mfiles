%% mean_curr_stress_crit

ncclear
clc
mp.urlC = ...
   'http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2009_da/his';
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';

nc=cfdataset(mp.urlC);
timeC = nc.time('ocean_time');
lon = nc.data('lon_rho');
lat = nc.data('lat_rho');

cin = find(timeC >= datenum(2010,04,29,00,00,00) & ...
    timeC <= datenum(2011,05,03,00,00,00));
timeC = timeC(cin);   %Cushion

var = nc.variable('ubar');
u = double(squeeze(var.data(cin,:,:)));
Ucur = u2rho_3d(u);
var = nc.variable('vbar');
v = double(squeeze(var.data(cin,:,:)));
Vcur = v2rho_3d(v);
var = nc.variable('angle');
angle_r = var.data(:);
clear u v var nc

nc=cfdataset(mp.urlH);
time = nc.time('time');
var = nc.variable('tauwc');
test = var.data(10,:,:); clear var
wPoints = find(~isnan(test)); clear test

mean_curr = NaN([length(time) size(lon,1) size(lon,2)]);

for pp = 1:length(wPoints)
    disp(['On pt = ' num2str(pp) ' of ' num2str(length(wPoints))]);
    ww = wPoints(pp);
    if all(isnan(Ucur(:,ww))) || all(isnan(Vcur(:,ww))), error('What?!?'), end
    u2 = smart_interp(timeC,Ucur(:,ww),time,3);
    v2 = smart_interp(timeC,Vcur(:,ww),time,3);
    Umag = sqrt(u2.^2 + v2.^2);
    phic = atan2(v2,u2) + angle_r(ww);
    mean_curr(:,ww) = (Umag.*cos(phic)) + 1i.*(Umag.*sin(phic));
end
clear phic Umag u2 v2 ww wPoints pp nc
clear Ucur Vcur angle_r cin test timeC
save depth_average_current