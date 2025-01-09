%% Uses actual grain size distribution

%1.  Use kN = 0.5 cm for the actual combined wave/current stress
%calculation
%2.  For each grid cell, determine if there is any sediment texture data;
%take the mean, standard deviation, and number of observations (of the
%mean)
%3.  Calculate a critical threshold value (save, spatially resolved)
%4.  Calculate the percentage of time above threshold


tic

%Pull grain size data
gFile = fullfile('C:\Users\sdalyander\Documents\SedimentTexture\FromBrad', ...
    'ecstdb2005.xls');
[nums,tex] = xlsread(gFile);
mc = nums(:,strcmp(tex(1,:),'MEAN'));
%Convert phi -> m
mc = (2.^(-1*mc))/1000;
my = nums(:,strcmp(tex(1,:),'LATITUDE'));
mx = nums(:,strcmp(tex(1,:),'LONGITUDE'));
p = [mx(:) my(:)]; clear mx my nums tex

%Stress observations
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN';

%% Open the dataset
nc = ncdataset(mp.urlH);
time = nc.time('time');
lon = nc.data('lon');
lat = nc.data('lat');
gtime = datevec(time);

%Pull down one time step as a mask
test = squeeze(double(nc.data('tauwc',[100 1 1],[100 size(lon,1) size(lon,2)])));

%Set up the empties
mean_gsize = NaN(size(lon));
std_gsize = NaN(size(lon));
num_gsize = NaN(size(lon));
crit_stress = NaN(size(lon));

rmpath('C:\Users\sdalyander\Documents\MATLAB\m_cmg\trunk\tri\')

%Borrowed from pcolorjw.m
[m n] = size(lon);
x = lon; y = lat;
x = [ x(:,1)  0.5*(x(:,1:n-1) + x(:,2:n))  x(:,n)];
y = [ y(:,1)  0.5*(y(:,1:n-1) + y(:,2:n))  y(:,n)];
x = [ x(1,:); 0.5*(x(1:m-1,:) + x(2:m,:)); x(m,:)];
y = [ y(1,:); 0.5*(y(1:m-1,:) + y(2:m,:)); y(m,:)];

%% Run through
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
            mean_gsize(ii,jj) = nanmean(mc(isIn));
            std_gsize(ii,jj) = nanstd(mc(isIn));
            [~,crit_stress(ii,jj)] = pmsoulsby(mean_gsize(ii,jj),1);
        end
    end
end


% save grid_sed_texture gFile mp mean_gsize std_gsize num_gsize crit_stress lon lat

%% Set up divisions
nn = 50;
jjL = unique([1:nn:size(test,1) size(test,1)]);
iiL = unique([1:nn:size(test,2) size(test,2)]);
clear nn

%stress types, analysis time periods
ttypes = {'wc'; 'w'; 'c'; 't'; 'r'};
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

mbed_mobil_act = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]); %season, grid

%% Run the loop
for jj = 1:length(jjL)-1
    for ii = 1:length(iiL)-1
        disp(['Data set ' num2str(jj) ',' num2str(ii) ' of ' ...
            num2str(length(jjL)-1) ',' num2str(length(iiL)-1)])
        
        jt = jjL(jj):jjL(jj+1)-1;
        it = iiL(ii):iiL(ii+1)-1;
        
        disp('Loading data')
        tauwc = nc.data('tauwc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        tauw = nc.data('tauw',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        tauc = nc.data('tauc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        taut = nc.data('tauc_tide',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        taur = nc.data('tauc_res',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        
        disp('Analyzing data')
        
        %Cycle through
        for mm = 1:length(mlists)
            inThis = find(ismember(gtime(:,2),mlists{mm}));
            for jj2 = 1:length(jt)
                for ii2 = 1:length(it)
                    %Set 5:  Geologic
                    if isnan(mean_gsize(jt(jj2),it(ii2)))
                        continue
                    end
                    for tt = 1:length(ttypes)
                        eval(['[mbed_mobil_act(tt,mm,jt(jj2),it(ii2))] = ' ...
                            'pmsoulsby(mean_gsize(jt(jj2),it(ii2)),tau' ...
                            ttypes{tt} '(inThis,jj2,ii2));'])
                    end
                end
            end
        end
    end
    clear tauwc tauw tauc tauw_wc tauc_wc taut_wc taur_wc taut taur
end
% close(nc), clear nc
clear nc
clear gtime test testCut iiL jjL ii2 jj2 ii jj tt mm Fs T nn it jt mp
clear nodes n p m inThis isIn wc_95 mean_wc x y mc
disp('Be sure to save results!!')
toc

