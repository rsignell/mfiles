%identify_crit_stress_from_obs

ncclear, clc
% load analysis1_gen_geo

sList = [10:10:90];

%Phi set up
phi = -5:1:11;
gsizel = (2.^(-1*phi))/1000;
phi = fliplr(phi); %smallest to largest (see below)
gsizel = fliplr(gsizel);

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

gFile2 = fullfile('C:\Users\sdalyander\Documents\SedimentTexture\FromBrad', ...
    'EN486_texture_usgs_bbutm.xls');
[nums,tex] = xlsread(gFile2);
first = find(strcmp(tex(1,:),'PHIm5'));
last = find(strcmp(tex(1,:),'PHI_11'));
mc2 = nums(:,first:last);
mc2 = fliplr(mc2);    %Go smallest to largest
mc2(isnan(mc2)) = 0;
mc2(mc2==-9999) = 0;
sum_mc2 = sum(mc2,2);
bads = find(sum_mc2 > 105 | sum_mc2 < 95);
mc2(bads,:) = NaN;
my2 = nums(:,strcmp(tex(1,:),'LATITUDE'));
mx2 = nums(:,strcmp(tex(1,:),'LONGITUDE'));
p2 = [mx2(:) my2(:)]; clear mx my nums tex
p = [p; p2]; clear p2
mc = [mc; mc2]; clear mc2
sum_mc = [sum_mc; sum_mc2]; clear sum_mc2 first last bads

% Open the dataset
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SABGOM_SWAN7_point/sabgom_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\SABGOM\stress_analysis';
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



grain_crit = NaN(size(gsizel));
for gg = 1:length(gsizel)
    [~,grain_crit(gg)] = pmsoulsby(gsizel(gg),1);
end

rmpath('C:\Users\sdalyander\Documents\MATLAB\m_cmg\trunk\tri\')
num_gsize = NaN(size(lon));
perc_crit_smaller_obs_all = NaN([size(mc,1) length(sList)]);
perc_crit_smaller_obs = NaN([length(sList) size(lon,1) size(lon,2)]);
model_stress_obs_pts = NaN([size(mc,1) length(sList)]);
box_obs_in = NaN([size(mc,1) 2]);
for ss = 1:length(sList)
    disp(['On ss = ' num2str(ss) ' of ' num2str(length(sList))])
    if sList(ss) == 50
        sVar = squeeze(mmed(1,1,:,:));
    else
        eval(['sVar = squeeze(m' num2str(sList(ss),'%02.0f') '(1,1,:,:));'])
    end
    
    
    for ii = 1:size(lon,1)
        disp(['On ii = ' num2str(ii) ' of ' num2str(size(lon,1))])
        for jj = 1:size(lon,2)
            
            nodes = [x(ii,jj) y(ii,jj); x(ii,jj+1) y(ii,jj+1); ...
                x(ii+1,jj+1) y(ii+1,jj+1); x(ii+1,jj) y(ii+1,jj)];
            [isIn,isOn] = inpoly(p, nodes);
            isIn = isIn | isOn; clear isOn
            
            if isempty(find(isIn, 1))
                num_gsize(ii,jj) = 0;
                continue
            end
            
            gall = mc(isIn,:);
            gall = nanmean(gall,1);
            if all(isnan(gall)), continue, end
            gall = cumsum(gall);
            
            hstress = sVar(ii,jj);
            if isnan(hstress), continue, end
            hin = find(grain_crit <= hstress,1,'last');
            pts = find(isIn);
            if isempty(hin)
                num_gsize(ii,jj) = length(find(isIn));
                perc_crit_smaller_obs(ss,ii,jj) = 0;
                for gg = 1:length(pts)
                    perc_crit_smaller_obs_all(pts(gg),ss) = 0;
                    model_stress_obs_pts(pts(gg),ss) = hstress;
                    box_obs_in(pts(gg),1) = ii;
                    box_obs_in(pts(gg),2) = jj;
                end
                continue
            end
            
            num_gsize(ii,jj) = length(find(~isnan(mc(isIn,1))));
            perc_crit_smaller_obs(ss,ii,jj) = gall(hin);
            for gg = 1:length(pts)
                box_obs_in(pts(gg),1) = ii;
                box_obs_in(pts(gg),2) = jj;
                model_stress_obs_pts(pts(gg),ss) = hstress;
                gall = mc(pts(gg),:);
                if all(isnan(gall)), continue, end
                gall(isnan(gall)) = 0;
                gall = cumsum(gall);
                perc_crit_smaller_obs_all(pts(gg),ss) = gall(hin);
                
            end
        end
    end
end

sList = 10:10:90;
gsize = NaN([size(p,1) length(sList)]);
cstress_grains = NaN([size(p,1) length(sList)]);
for ss = 1:length(sList)
    for gg = 1:size(mc,1)
        thisg = mc(gg,:);
        thisg(isnan(thisg)) = 0;
        if all(thisg == 0), clear thisg, continue, end
        csum = cumsum(thisg); clear thisg
        fover = find(csum >= sList(ss),1,'first'); clear csum
        gsize(gg,ss) =  gsizel(fover);
        [~,cstress_grains(gg,ss)] = pmsoulsby(gsize(gg,ss),1);
    end
end



clear first gall gg hin hstress ii jj isIn last m m05 m10 m25 m75 m95 gsize
clear mbed_mobil mlists mmean mmed mspread mstd n nodes sVar ss
clear tauc_perc taur_perc taut_perc tauw_perc time ttypes x y pts
clear m20 m30 m40 m60 m70 m80 m90 mmed
save perc_grainsize_crit_smaller_percentiles

% 
% %% Plot
% edges = [0:5:100 1000];
% edges2 = [0:5:100];
% close all
% figure
% for ss = 1:length(sList)
%     %     sVar = squeeze(perc_below(ss,:,:));
%     sVar = perc_crit_smaller_obs_all(:,ss);
%     subplot(3,3,ss)
%     n = histc(sVar(:),edges);
%     n(end-2) = n(end-2)+n(end-1);
%     n(end-1) = 0;
%     n(end) = [];
%     bar(edges2,n,'histc')
%     ylim([0 1000])
%     xlim([0 100])
%     title({'Histogram of Percent Sediment'; ['\tau_c_r_i_t < ' num2str(sList(ss)) 'th Percentile Stress']})
% end
