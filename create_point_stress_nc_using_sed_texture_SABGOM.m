%% create_point_stress_nc_using_sed_texture

% cd('C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\point_stress')
ncclear, clc

nStart = 1;

%% Phi class set up
phi = -5:1:11;
gsizel = (2.^(-1*phi))/1000;
phi = fliplr(phi); %smallest to largest (see below)
gsizel = fliplr(gsizel);

%% Pull grain size data, file 1
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
locs = [mx(:) my(:)]; clear mx my nums tex

%% Pull grain size data, file 2
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
locs = [locs; p2]; clear p2
mc = [mc; mc2]; clear mc2
sum_mc = [sum_mc; sum_mc2]; clear sum_mc2 first last bads
csum = cumsum(mc,2);
clear mx2 my2

%%  Data for the files
% mp.history = ['Grain size data from ' gFile ' and ' gFile2 '; Waves from '...
%     'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CompGrid/swan_agg.ncml; ' ...
%     'Currents from Espresso forecast as of 11/30/2011'];
mp.history = ['Grain size data from ' gFile ' and ' gFile2 '; Waves from '...
    'SWAN 7 Grid, Currents from SABGOM forecast'];
clear gFile gFile2
mp.type = 'Stress time series using point measurements of sediment texture for physical roughness.';
mp.NL = size(mc,1); %Number of sediment observations
mp.CL = size(mc,2); %Number of phi classes

mp.swan = 'C:\Users\sdalyander\Documents\Waves\SWAN_7grid\SWAN_hind_kNpt05\';
cd(mp.swan)

%% Create all the empty files, one per month
wList = [1005:1012 1101:1104];
for ii = nStart:length(wList)
    mp.wrun = wList(ii);
    data = load(fullfile(mp.swan,['wrun' num2str(mp.wrun) '_sabgom_ubot.mat']),'time');
    timeW = data.time; clear data
    mp.tL = length(timeW);
    swrun = num2str(mp.wrun);
    yyyy = ['20' swrun(1:2)];
    mm = swrun(3:4); clear swrun
    fileN{ii} = ['sab_point_' yyyy mm '.nc'];
    status = create_point_tau_output(fileN{ii},mp);
    
    myTime =  timeW - datenum(1858,11,17);
    
    nc = netcdf(fileN{ii},'write');
    nc{'time'}(:) = myTime; clear myTime
    nc{'lon'}(:) = NaN(mp.NL,1);
    nc{'lat'}(:) = NaN(mp.NL,1);
    nc{'kN'}(:) = NaN(mp.NL,1);
    nc{'tau_crit'}(:) = NaN(mp.NL,1);
    nc{'zr'}(:,:) = NaN(mp.tL,mp.NL);
    nc{'tauwc'}(:,:) = NaN(mp.tL,mp.NL);
    close(nc)
end
clear status ii nc

%% ROMS stuff, pull down
%Current parameters
mp.urlG = ...
    'http://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SABGOM_2011/sabgom_grd_H.nc';
mp.urlCU = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SABGOM_2011/sabgom_u_agg.ncml';
mp.urlCV = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SABGOM_2011/sabgom_v_agg.ncml';
mp.urlZ = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SABGOM_2011/sabgom_zeta_agg.ncml';
ncG = mDataset(mp.urlG);
ncU = mDataset(mp.urlCU);
lonC = ncG{'lon_rho'}(:);
latC = ncG{'lat_rho'}(:);
maskC = ncG{'mask_rho'}(:);
wPoints = find(maskC == 1);
timeC_all = nj_time(ncU, 'u');
angle_r = ncG{'angle'}(:);
[mp.RL,mp.CL] = size(lonC);

%Pull parameters for depth calculation
%Follows roms_get_grid
Cs_r =  [-0.9593; -0.8858; -0.8213; -0.7646; -0.7145; -0.6700; -0.6300; ...
    -0.5936; -0.5599; -0.5280; -0.4973; -0.4668; -0.4360; -0.4043; -0.3713; ...
    -0.3372; -0.3020; -0.2665; -0.2315; -0.1980; -0.1666; -0.1383; -0.1133; ...
    -0.0918; -0.0737; -0.0586; -0.0463; -0.0362; -0.0281; -0.0216; -0.0163; ...
    -0.0120; -0.0085; -0.0055; -0.0031; -0.0010]; %Taken from SABGOM THREDDS
sc_r = ncU{'s_rho'}(:);
h = ncG{'h'}(:);    %Depth as positive
hc = 5; %SABGOM THREDDS
scmCshc = (sc_r-Cs_r)*hc;
z_rB = repmat(scmCshc,[1 length(h(:))]) + Cs_r*h(:)';
close(ncG), close(ncU)
clear mVar hc N ncG ncU ncZ

% %% Set up the boxes for point location on grid & find the points
% %for SABGOM, already have SWAN interpolated to current grid
% [m n] = size(lonC);
% xC = lonC; yC = latC;
% xC = [ xC(:,1)  0.5*(xC(:,1:n-1) + xC(:,2:n))  xC(:,n)];
% yC = [ yC(:,1)  0.5*(yC(:,1:n-1) + yC(:,2:n))  yC(:,n)];
% xC = [ xC(1,:); 0.5*(xC(1:m-1,:) + xC(2:m,:)); xC(m,:)];
% yC = [ yC(1,:); 0.5*(yC(1:m-1,:) + yC(2:m,:)); yC(m,:)];
% clear m n
% rmpath('C:\Users\sdalyander\Documents\MATLAB\m_cmg\trunk\tri\')
% Cboxes = NaN(size(locs));
% isGood = find(maskC(:) == 1);
% for bb = 1:length(isGood)
%     [ii,jj] = ind2sub(size(lonC),isGood(bb));
%     nodes = [xC(ii,jj) yC(ii,jj); xC(ii,jj+1) yC(ii,jj+1); ...
%         xC(ii+1,jj+1) yC(ii+1,jj+1); xC(ii+1,jj) yC(ii+1,jj)];
%     [isIn,isOn] = inpoly(locs, nodes);
%     isIn = isIn | isOn; clear isOn
%     Cboxes(isIn,1) = ii;
%     Cboxes(isIn,2) = jj;
% end
%
load(fullfile('C:\Users\sdalyander\Documents\StressAnalysis\SABGOM', ...
    'point_boxes'))
Cinds = sub2ind(size(h),Cboxes(:,1),Cboxes(:,2));

%% Other needed things
grain_crit = NaN(size(gsizel));
for gg = 1:length(gsizel)
    [~,grain_crit(gg)] = pmsoulsby(gsizel(gg),1);
end
clear gg
cohesive_split = find(gsizel < 0.004e-3,1,'last');
phi_cohesive = (1/16)*10^-3;    %4 Phi, split of sand and silt
[~,crit_cohesive] = pmsoulsby(phi_cohesive,1);


%% Go sample by sample, , and fill the files
zrMax = 5;
done = zeros(mp.NL,1);

Cboxes(Cboxes(:,1) == size(h,1),:) = NaN;
Cboxes(Cboxes(:,2) == size(h,2),:) = NaN;

Cinds = sub2ind(size(lonC),Cboxes(:,1),Cboxes(:,2));

%%
for ww = nStart:length(wList)    %Loop through the months
    disp(['On ww = ' num2str(ww) ' of ' num2str(length(wList))])
    
    %Pull the data
    disp('Pulling in the wave model output')
    mp.wrun = wList(ww);
    ubData = load(fullfile(mp.swan,['wrun' num2str(mp.wrun) '_sabgom_ubot.mat']),'ubot','time');
    tbData = load(fullfile(mp.swan,['wrun' num2str(mp.wrun) '_sabgom_tmbot.mat']),'tmbot');
    wdData = load(fullfile(mp.swan,['wrun' num2str(mp.wrun) '_sabgom_bdir.mat']),'bdir');
    
    timeW = ubData.time;
    ubr = ubData.ubot;
    tbr = tbData.tmbot;
    wdir = wdData.bdir; %10/02/2011 change to dirbot
    clear ubData tbData wdData
    
    isIn = find((timeC_all >= min(timeW) - 0.5/24) & (timeC_all <= max(timeW) + 0.5/24));
    timeC = timeC_all(isIn);
    
    [timeC,II] = unique(timeC);
    
    for nn2 = 1:mp.NL
        if done(nn2) || isnan(Cboxes(nn2,1)), continue, end
        
        disp(['Analyzing nn2 = ' num2str(nn2) ' of ' num2str(mp.NL)])
        
        %Pull the wave data
        ubrA = ubr(:,Cboxes(nn2,1),Cboxes(nn2,2));
        tbrA = tbr(:,Cboxes(nn2,1),Cboxes(nn2,2));
        wdirA = wdir(:,Cboxes(nn2,1),Cboxes(nn2,2));
        wdirA = (wdirA*pi/180)*-1 + (3*pi/2); %Degrees relative to North, direction waves coming from to radians rel. East, direction to
        
        if all(isnan(ubrA)), continue, end
        
        %Water level
        ncZ = mDataset(mp.urlZ);
        zeta = double(ncZ{'zeta'}(isIn,Cboxes(nn2,1),Cboxes(nn2,2)));
        close(ncZ), clear ncZ
        zeta = zeta(II);
        
%         disp('Extracting depths')
        zrO = NaN(1,length(timeW));
        for ii = 1:length(timeW)
            zeta1 = zeta(ii);
            z_r = z_rB(1,Cinds(nn2)) + scmCshc(1)*[zeta1./h(Cinds(nn2))]' + (1+Cs_r(1))*zeta1';
            zrO(ii) = z_r + h(Cinds(nn2));
        end
        
        ll = 1;
        ncCU = mDataset(mp.urlCU);
        ncCV = mDataset(mp.urlCV);
        u = double(squeeze(ncCU{'u'}(isIn,ll,Cboxes(nn2,1),Cboxes(nn2,2)-1:Cboxes(nn2,2))));
        v = double(squeeze(ncCV{'v'}(isIn,ll,Cboxes(nn2,1)-1:Cboxes(nn2,1),Cboxes(nn2,2))));
        close(ncCU), close(ncCV), clear ncCU ncCV
        curr = nanmean(u,2) + 1i.*nanmean(v,2); clear u v
        
        if all(isnan(curr)), continue, end
        
        isThis = find(Cinds == Cinds(nn2));
        
        for nn3 = 1:length(isThis)
            nn = isThis(nn3);
            
            zrN = NaN(length(timeW),1);
            tauwc = NaN(length(timeW),1);
            kN = NaN;
            tau_crit = NaN;
            
            if isnan(csum(nn,end)),continue,end
            
            %Check for cohesive
            if csum(nn,cohesive_split) >=7.5
                kN = phi_cohesive;
                tau_crit = crit_cohesive;
            else
                kN = gsizel(find(csum(nn,:) >= 50,1));
                tau_crit = grain_crit(find(csum(nn,:) >= 50,1));
            end
            
%             disp('Starting time analysis')
            for tt = 1:length(timeW)
                %Check current model depth
                zr = zrO(tt);
                if zr > zrMax, continue, end
                
                Umag = sqrt(real(curr(tt)).^2 + imag(curr(tt)).^2 + eps);
                phic = atan2(imag(curr(tt)),real(curr(tt))) + angle_r(Cboxes(nn,1),Cboxes(nn,2));
                phiwc = abs(wdirA(tt) - phic);
                
                
                m = m94(ubrA(tt),2*pi/tbrA(tt),Umag,zr,phiwc,kN,0);
                
                while m.dwc > zr
                    ll = ll+1;
                    ncCU = mDataset(mp.urlCU);
                    ncCV = mDataset(mp.urlCV);
                    gg = value2Index(timeC_orig,timeW(tt),0.5);
                    u = double(squeeze(ncCU{'u'}(gg,ll,:,:)));
                    Ucur = u2rho_2d(u);
                    v = double(squeeze(ncCV{'v'}(gg,ll,:,:)));
                    Vcur = v2rho_2d(v);
                    ii_w = Cboxes(nn,1); jj_w = CBoxes(nn,2);
                    Umag = sqrt(Ucur(ii_w,jj_w).^2 + Vcur(ii_w,jj_w).^2 + eps);
                    phic = atan2(Vcur(ii_w,jj_w),Ucur(ii_w,jj_w)) + angle_r(ii_w,jj_w); %Rotate
                    
                    zeta1 = zeta(tt);
                    z_r = z_rB(ll,Cinds(nn2)) + scmCshc(ll)*[zeta1./h(Cinds(nn2))]' + (1+Cs_r(ll))*zeta1';
                    zr = z_r + h(Cboxes(nn,1),Cboxes(nn2));
                    
                    phiwc(ww) = abs(wdirA(tt) - phic);
                    m = m94(ubrA(tt),2*pi/tbrA(tt),Umag(tt),zr,phiwc,kN,0);
                    clear u v curr_out zeta1
                    close(ncCU), close(ncCV), clear ncCU ncCV
                end %Ends while loop checking zr
                zrN(tt) = zr;
                tauwc(tt) = 1025*m.ustrr.^2;
            end %Ends time series loop
            
            ncR = netcdf(fileN{ww},'write');
            ncR{'lon'}(nn) = locs(nn,1);
            ncR{'lat'}(nn) = locs(nn,2);
            ncR{'kN'}(nn) = kN;
            ncR{'tau_crit'}(nn) = tau_crit;
            ncR{'zr'}(:,nn) = zrN;
            ncR{'tauwc'}(:,nn) = tauwc;
            close(ncR)
            done(nn) = 1;
        end %Ends loop of points with same boxes
    end %Ends box loop
end %Ends month loop

