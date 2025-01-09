%% spatial_stress_analysis_p3

%Process-based analysis
%CAN'T RUN ON MORE THAN A YEAR OF DATA for seasonal analysis
ncclear, clc
tic
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SABGOM_SWAN7/sabgom_agg.ncml';
mp.mdir = 'C:\Users\sdalyander\Documents\StressAnalysis\SABGOM\stress_analysis';

%% Open the dataset
% nc = mDataset(mp.urlH);
% time = nj_time(nc,'tauwc');
% lon = double(nc{'lon'}(:));
% lat = double(nc{'lat'}(:));
nc = ncdataset(mp.urlH);
time = nc.time('time');
lon = nc.data('lon');
lat = nc.data('lat');

gtime = datevec(time);

%Pull down one time step as a mask
% test = nc{'tauwc'}(100,:,:);
test = squeeze(double(nc.data('tauwc',[100 1 1],[100 size(lon,1) size(lon,2)])));

%Set up divisions
nn = 50;
jjL = unique([1:nn:size(test,1) size(test,1)]);
iiL = unique([1:nn:size(test,2) size(test,2)]);
clear nn


load(fullfile('C:\Users\sdalyander\Documents\Winds_Models', 'wind_buoys'))

%stress types, analysis time periods
ttypes = {'wc'; 'w'; 'c'; 't'; 'r'};
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

%Set 4:Processed based
xcor_wind_stress = NaN([length(buoyList) length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
xcor_wind_u = NaN([length(buoyList) length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
xcor_wind_v = NaN([length(buoyList) length(ttypes) length(mlists) size(lon,1) size(lon,2)]);

wind_tau2 = NaN([length(time) length(buoyList)]);
for bb = 1:length(buoyList)
    wind_tau2(:,bb) = smart_interp(timeNAM,real(wind_tau(:,bb,2)),time,5) + ...
        1i*smart_interp(timeNAM,imag(wind_tau(:,bb,2)),time,5);
end
    

%% Run the loop
for jj = 1:length(jjL)-1
    for ii = 1:length(iiL)-1
        disp(['Data set ' num2str(jj) ',' num2str(ii) ' of ' ...
            num2str(length(jjL)-1) ',' num2str(length(iiL)-1)])
        
        jt = jjL(jj):jjL(jj+1)-1;
        it = iiL(ii):iiL(ii+1)-1;
        
        disp('Loading data')
        %         tauwc = squeeze(double(nc{'tauwc'}(:,jt,it)));
        %         tauw = squeeze(double(nc{'tauw'}(:,jt,it)));
        %         tauc = squeeze(double(nc{'tauc'}(:,jt,it)));
        tauwc = nc.data('tauwc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        tauw = nc.data('tauw',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        tauc = nc.data('tauc',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        taut = nc.data('tauc_tide',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        taur = nc.data('tauc_res',[1 min(jt) min(it)],[length(time) max(jt) max(it)]);
        
        disp('Analyzing data')
        
        testCut = test(jt,it);
        
        %Cycle through
        for jj2 = 1:length(jt)
            for ii2 = 1:length(it)
                if isnan(testCut(jj2,ii2)), continue, end
                
                for mm = 1:length(mlists)
                    inThis = find(ismember(gtime(:,2),mlists{mm}));
                    
                    tauwc_it = tauwc(inThis,jj2,ii2);
                    tauw_it = tauw(inThis,jj2,ii2);
                    tauc_it = tauc(inThis,jj2,ii2);
                    taut_it = taut(inThis,jj2,ii2);
                    taur_it = taur(inThis,jj2,ii2);
                    
                    
                    for tt = 1:length(ttypes)
                        
                        for bb = 1:length(buoyList)
                            eval(['ttau = tau' ttypes{tt} '_it;']);
                            bstress = wind_tau2(inThis,bb);
                            [R,P] = corrcoef(abs(bstress),ttau,'rows','pairwise');
                            if P(1,2) > 0.05    %Insignificant correlation
                                xcor_wind_stress(bb,tt,mm,jt(jj2),it(ii2)) = 0;
                            else    %Significant correlation
                                xcor_wind_stress(bb,tt,mm,jt(jj2),it(ii2)) = R(1,2);
                            end
                            
                                                        [R,P] = corrcoef(real(bstress),ttau,'rows','pairwise');
                            if P(1,2) > 0.05    %Insignificant correlation
                                xcor_wind_u(bb,tt,mm,jt(jj2),it(ii2)) = 0;
                            else    %Significant correlation
                                xcor_wind_u(bb,tt,mm,jt(jj2),it(ii2)) = R(1,2);
                            end
                            
                                                        [R,P] = corrcoef(imag(bstress),ttau,'rows','pairwise');
                            if P(1,2) > 0.05    %Insignificant correlation
                                xcor_wind_v(bb,tt,mm,jt(jj2),it(ii2)) = 0;
                            else    %Significant correlation
                                xcor_wind_v(bb,tt,mm,jt(jj2),it(ii2)) = R(1,2);
                            end
                        end
                    end
                end
            end
        end
        clear tauwc tauw tauc tauw_wc tauc_wc taut_wc taur_wc taut taur
        clear tauwc_it tauw_it tauc_it taut_it taur_it inThis
    end
end
clear nc
clear gtime test testCut iiL jjL ii2 jj2 ii jj tt mm Fs T nn it jt mp
clear ind iw jw
disp('Be sure to save results!!')
toc
clear P R ans bb blon bhull c data_source lags loc2 u ttau
