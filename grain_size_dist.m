% grain_size_dist


ncclear, clc, close all
%Pull grain size data
%Pull grain size data
gFile = fullfile('C:\Users\sdalyander\Documents\SedimentTexture\FromLarry', ...
    'ecstdb2011.xls');
[nums,tex] = xlsread(gFile);
first = find(strcmp(tex(1,:),'PHIM5'));
last = find(strcmp(tex(1,:),'PHI_11'));
mc = nums(:,first:last);
mc = fliplr(mc);    %Go smallest to largest
mc(isnan(mc)) = 0;
mc(mc==-9999) = 0;
sum_mc = sum(mc,2);
bads = find(sum_mc > 105 | sum_mc < 95);
mc(bads,:) = NaN;
my = nums(:,strcmp(tex(1,:),'LATITUDE'));
mx = nums(:,strcmp(tex(1,:),'LONGITUDE'));
p = [mx(:) my(:)]; clear mx my nums tex

thresh = 10;    %percent

gsize = NaN([size(p,1) 1]);
for gg = 1:size(mc,1)
    thisg = mc(gg,:);
    thisg(isnan(thisg)) = 0;
    if all(thisg == 0), clear thisg, continue, end
    csum = cumsum(thisg); clear thisg
    fover = find(csum >= thresh,1,'first'); clear csum 
    gsize(gg) =  gsizel(fover);
end

%Stress observations
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';

% Open the dataset
nc = ncdataset(mp.urlH);
lon = nc.data('lon');
lat = nc.data('lat');
clear nc

%Borrowed from pcolorjw.m
[m n] = size(lon);
x = lon; y = lat;
x = [ x(:,1)  0.5*(x(:,1:n-1) + x(:,2:n))  x(:,n)];
y = [ y(:,1)  0.5*(y(:,1:n-1) + y(:,2:n))  y(:,n)];
x = [ x(1,:); 0.5*(x(1:m-1,:) + x(2:m,:)); x(m,:)];
y = [ y(1,:); 0.5*(y(1:m-1,:) + y(2:m,:)); y(m,:)];

%% Run through
%Set up the empties
mean_gsize = NaN(size(lon));
std_gsize = NaN(size(lon));
num_gsize = NaN(size(lon));
crit_stress = NaN(size(lon));

for ii = 1:size(lon,1)
    disp(['On ii = ' num2str(ii) ' of ' num2str(size(lon,1))])
    for jj = 1:size(lon,2)
        
        nodes = [x(ii,jj) y(ii,jj); x(ii,jj+1) y(ii,jj+1); ...
            x(ii+1,jj+1) y(ii+1,jj+1); x(ii+1,jj) y(ii+1,jj)];
        [isIn,isOn] = inpoly(p, nodes);
        isIn = isIn | isOn; clear isOn
        
        if isempty(isIn)
            num_gsize(ii,jj) = 0;
        else
            num_gsize(ii,jj) = length(find(isIn));
            mean_gsize(ii,jj) = nanmean(gsize(isIn));
            std_gsize(ii,jj) = nanstd(gsize(isIn));
            [~,crit_stress(ii,jj)] = pmsoulsby(mean_gsize(ii,jj),1);
        end
    end
end
eval(['mean_' num2str(thresh) 'p_gsize = mean_gsize;'])
eval(['std_' num2str(thresh) 'p_gsize = std_gsize;'])
eval(['crit_' num2str(thresh) 'p_stress = crit_stress;'])
clear mean_gsize std_gsize crit_stress ii jj isIn x y nodes first fover gg
clear last m n 
% disp('Don''t forget to save results!!')
eval(['save grid_sed_texture_' num2str(thresh) 'p_thresh'])
