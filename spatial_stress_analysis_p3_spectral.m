%% spatial_stress_analysis_p3

%Process-based analysis
%CAN'T RUN ON MORE THAN A YEAR OF DATA for seasonal analysis
ncclear
clc
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



%stress types, analysis time periods
ttypes = {'wc'; 'w'; 'c'; 't'; 'r'};
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

%Set 4:Processed based

spec_analysis_high = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);
spec_analysis_low = NaN([length(ttypes) length(mlists) size(lon,1) size(lon,2)]);

%time cut off for spectral analysis
Tcut = 33;  %hours
% Hs = spectrum.periodogram;
Hs = spectrum.welch;
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
                        eval(['ttau = tau' ttypes{tt} '_it;']);
                         ttau = ttau - nanmean(ttau);
                         try
                        hmss = msspectrum(Hs,ttau,'Fs',1);
                         catch
                             continue
                         end
                        f = hmss.Frequencies;
                        pow = hmss.Data;
                        low = (f<=1/Tcut);
                        high = (f>1/Tcut);
                        highfrac = trapz(f(high),pow(high));
                        lowfrac = trapz(f(low),pow(low));
                        spec_analysis_high(tt,mm,jt(jj2),it(ii2)) = highfrac;
                        spec_analysis_low(tt,mm,jt(jj2),it(ii2)) = lowfrac;
                    end
                    clear frac high low pow f Y ttau R P L t nfft
                end
            end
        end
        clear tauwc tauw tauc tauw_wc tauc_wc taut_wc taur_wc taut taur
        clear tauwc_it tauw_it tauc_it taut_it taur_it inThis
    end
end

clear nc nc2 hmss Hs 
clear gtime test testCut iiL jjL ii2 jj2 ii jj tt mm Fs T nn it jt mp
clear ind iw jw
disp('Be sure to save results!!')
toc

