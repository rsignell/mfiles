%% spatial_stress_analysis_p3

%%IS NOT YET OPERATIONAL!  NEEDS TO BE CORRECTED FOR LON/LAT OF WIND GRID

%Process-based analysis
%CAN'T RUN ON MORE THAN A YEAR OF DATA for seasonal analysis
ncclear
tic
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CompGrid/swan_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\WavesOnly';

%% Open the dataset
% nc = mDataset(mp.urlH);
% time = nj_time(nc,'tauwc');
% lon = double(nc{'lon'}(:));
% lat = double(nc{'lat'}(:));
nc = ncdataset(fullfile('C:\Users\sdalyander\Documents\COAWST\EC30Day\grids', ...
    'USeast_grd17_psd3.nc'));
lon = nc.data('lon_rho');
lat = nc.data('lat_rho');
mask = nc.data('mask_rho');
clear nc

nc = ncdataset(mp.urlH);
time = nc.time('time');
gtime = datevec(time);

%Pull down one time step as a mask
% test = nc{'tauwc'}(100,:,:);
test = squeeze(double(nc.data('tauw',[100 1 1],[100 size(lon,1) size(lon,2)])));

%Set up divisions
nn = 50;
jjL = unique([1:nn:size(test,1) size(test,1)]);
iiL = unique([1:nn:size(test,2) size(test,2)]);
clear nn

%Winds
loc2 = 'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/forcings_winds_USEast/wind_agg.ncml';
ncWind=mDataset(loc2);
timeWind = nj_time(ncWind,'Uwind');
sWind = value2Index(timeWind,min(time));
eWind = value2Index(timeWind,max(time));
timeWind = timeWind(sWind:eWind);
lonWind = 

%stress types, analysis time periods
ttypes = {'wc'; 'w'; 'c'; 't'; 'r'};
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

%Set 4:Processed based
% xcor_wave_curr = NaN([length(mlists) size(lon,1) size(lon,2)]); %Cross correlation waves x currents
% xcor_wave_res = NaN([length(mlists) size(lon,1) size(lon,2)]);  %Cross correlation waves x currents/no tides
% xcor_wind_all = NaN([length(mlists) size(lon,1) size(lon,2)]);
xcor_wind_wave = NaN([length(mlists) size(lon,1) size(lon,2)]);
% xcor_wind_curr = NaN([length(mlists) size(lon,1) size(lon,2)]);
% xcor_wind_tide = NaN([length(mlists) size(lon,1) size(lon,2)]);
% xcor_wind_res = NaN([length(mlists) size(lon,1) size(lon,2)]);


%% Run the loop
for jj = 1:length(jjL)-1
    for ii = 1:length(iiL)-1
        disp(['Data set ' num2str(jj) ',' num2str(ii) ' of ' ...
            num2str(length(jjL)-1) ',' num2str(length(iiL)-1)])
        
        jt = jjL(jj):jjL(jj+1)-1;
        it = iiL(ii):iiL(ii+1)-1;
        if all(all(mask(jt,it)==0)),continue,end
        
        disp('Loading data')
        %         tauwc = squeeze(double(nc{'tauwc'}(:,jt,it)));
        %         tauw = squeeze(double(nc{'tauw'}(:,jt,it)));
        %         tauc = squeeze(double(nc{'tauc'}(:,jt,it)));
        %         tauwc = nc.data('tauwc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        tauw = nc.data('tauw',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        if all(isnan(tauw(:))) | all(tauw(:) == 0), continue, end
%         tauc = nc.data('tauc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        %         taut = nc.data('tauc_tide',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
%         taur = nc.data('tauc_res',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        
        
        %Pull wind
        disp('Pulling wind')
        
        u_wind_all = ncWind{'Uwind'}(sWind:eWind,jt,it);
        v_wind_all = ncWind{'Vwind'}(sWind:eWind,jt,it);
        disp('Wind pulled, running analysis')
        
        disp('Analyzing data')
        
        testCut = test(jt,it);
        
        %Cycle through
        for jj2 = 1:length(jt)
            for ii2 = 1:length(it)
                if isnan(testCut(jj2,ii2)), continue, end
%                 disp('Interp wind')
                [u_tau,v_tau] = wstress(squeeze(u_wind_all(:,jj2,ii2)),...
                    squeeze(v_wind_all(:,jj2,ii2)));
                u_tau = smart_interp(timeWind,u_tau,time,5);
                v_tau = smart_interp(timeWind,v_tau,time,5);
                ws = sqrt(u_tau.^2 + v_tau.^2);
                clear u_tau v_tau u_wind v_wind jw iw ind
                for mm = 1%:length(mlists)
                    inThis = find(ismember(gtime(:,2),mlists{mm}));
                    
                    %                     tauwc_it = tauwc(inThis,jj2,ii2);
                    tauw_it = tauw(inThis,jj2,ii2);
%                     tauc_it = tauc(inThis,jj2,ii2);
                    %                     taut_it = taut(inThis,jj2,ii2);
%                     taur_it = taur(inThis,jj2,ii2);
                    
                    %
                    %                     [R,P] = corrcoef(tauw_it,tauc_it);
                    %                     if P(1,2) > 0.05    %Insignificant correlation
                    %                         xcor_wave_curr(mm,jt(jj2),it(ii2)) = 0;
                    %                     else    %Significant correlation
                    %                         xcor_wave_curr(mm,jt(jj2),it(ii2)) = R(1,2);
                    %                     end
                    %
                    %                     [R,P] = corrcoef(tauw_it,taur_it);
                    %                     if P(1,2) > 0.05    %Insignificant correlation
                    %                         xcor_wave_res(mm,jt(jj2),it(ii2)) = 0;
                    %                     else    %Significant correlation
                    %                         xcor_wave_res(mm,jt(jj2),it(ii2)) = R(1,2);
                    %                     end
                    %
                    %                     [R,P] = corrcoef(ws(inThis),tauwc_it);
                    %                     if P(1,2) > 0.05
                    %                         xcor_wind_all(mm,jt(jj2),it(ii2)) = 0;
                    %                     else
                    %                         xcor_wind_all(mm,jt(jj2),it(ii2)) = R(1,2);
                    %                     end
                    
                    [R,P] = corrcoef(ws(inThis),tauw_it);
                    if P(1,2) > 0.05
                        xcor_wind_wave(mm,jt(jj2),it(ii2)) = 0;
                    else
                        xcor_wind_wave(mm,jt(jj2),it(ii2)) = R(1,2);
                    end
                    
%                     [R,P] = corrcoef(ws(inThis),tauc_it);
%                     if P(1,2) > 0.05
%                         xcor_wind_curr(mm,jt(jj2),it(ii2)) = 0;
%                     else
%                         xcor_wind_curr(mm,jt(jj2),it(ii2)) = R(1,2);
%                     end
%                     %
%                     %                     [R,P] = corrcoef(ws(inThis),taut_it);
%                     %                     if P(1,2) > 0.05
%                     %                         xcor_wind_tide(mm,jt(jj2),it(ii2)) = 0;
%                     %                     else
%                     %                         xcor_wind_tide(mm,jt(jj2),it(ii2)) = R(1,2);
%                     %                     end
%                     
%                     [R,P] = corrcoef(ws(inThis),taur_it);
%                     if P(1,2) > 0.05
%                         xcor_wind_res(mm,jt(jj2),it(ii2)) = 0;
%                     else
%                         xcor_wind_res(mm,jt(jj2),it(ii2)) = R(1,2);
%                     end
                    
                end
            end
        end
        clear tauwc tauw tauc tauw_wc tauc_wc taut_wc taur_wc taut taur
        clear tauwc_it tauw_it tauc_it taut_it taur_it inThis
    end
end
% close(nc), clear nc
close(ncWind), clear ncNAM
clear nc
clear gtime test testCut iiL jjL ii2 jj2 ii jj tt mm Fs T nn it jt mp
 clear lonWind latWind sWind eWind P R ii_w jj_w ws timeWind u_wind_all v_wind_all
clear ind iw jw
disp('Be sure to save results!!')
toc

