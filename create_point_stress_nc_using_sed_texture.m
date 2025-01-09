%% create_point_stress_nc_using_sed_texture

cd('C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\point_stress')
ncclear, clc

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

%% Pull the SWAN metadata
nc = mDataset(fullfile('C:\Users\sdalyander\Documents\COAWST\EC30Day\grids', ...
    'USeast_grd17_psd3.nc'));
% maskSWAN = nc{'mask_rho'}(:);
lonSWAN = nc{'lon_rho'}(:);
% latSWAN = nc{'lat_rho'}(:);
% maskSWAN(~maskSWAN) = NaN;
close(nc)
nc = ncdataset('dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CompGrid/swan_agg.ncml');
timeSWAN = nc.time('time');
clear nc

%%  Data for the files
mp.history = ['Grain size data from ' gFile ' and ' gFile2 '; Waves from '...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CompGrid/swan_agg.ncml; ' ...
    'Currents from Espresso forecast as of 11/30/2011'];
clear gFile gFile2
mp.type = 'Stress time series using point measurements of sediment texture for physical roughness.';
mp.NL = size(mc,1); %Number of sediment observations
mp.CL = size(mc,2); %Number of phi classes


%% Create all the empty files, one per month
yyList = [repmat(2010,1,8) repmat(2011,1,4)];
mmList = [5:12 1:4];
gTime = datevec(timeSWAN);
for ii = 1:length(mmList)
    mp.tL = length(find(gTime(:,1) == yyList(ii) & gTime(:,2) == mmList(ii)));
    fileN{ii} = ['mab_point_' num2str(yyList(ii)) num2str(mmList(ii),'%02.0f') '.nc'];
    status = create_point_tau_output(fileN{ii},mp);
    
    myTime =  timeSWAN(gTime(:,1) == yyList(ii) & gTime(:,2) == mmList(ii)) - datenum(1858,11,17);
    
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
mp.urlC = ...
    'http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2009_da/his';
ncC = mDataset(mp.urlC);
% lonC = ncC{'lon_rho'}(:);
% latC = ncC{'lat_rho'}(:);
% maskC = ncC{'mask_rho'}(:);
timeC = nj_time(ncC, 'u');
angle_r = ncC{'angle'}(:);
h = ncC{'h'}(:);
Cs_r = ncC{'Cs_r'}(:);
Cs_w = ncC{'Cs_w'}(:);
sc_r = ncC{'s_rho'}(:);
sc_w = ncC{'s_w'}(:);
N = length(sc_r);
mVar = getVar(ncC, 'hc');
hc = getData(mVar);
if isempty(hc)
    disp('Assuming min(h) = hc');
    hc = min(min(h));
end
scmCshc = (sc_r-Cs_r)*hc;
z_rB = repmat(scmCshc,[1 length(h(:))]) + Cs_r*h(:)';
scmCshc_w = (sc_w-Cs_w)*hc;
z_wB = repmat(scmCshc_w,[1 length(h(:))]) + Cs_w*h(:)';
close(ncC)
clear hc N mVar sc_r sc_w

%% Pull down the information
load(fullfile('C:\Users\sdalyander\Documents\StressAnalysis\MAB_Espresso_SWAN\point_stress', ...
    'containing_boxes'))
SWANinds = sub2ind(size(lonSWAN),SWANboxes(:,1),SWANboxes(:,2));
Cinds = sub2ind(size(h),Cboxes(:,1),Cboxes(:,2));

% %% Set up the boxes for point location on grid & find the points
% [m n] = size(lonSWAN);
% xSWAN = lonSWAN; ySWAN = latSWAN;
% xSWAN = [ xSWAN(:,1)  0.5*(xSWAN(:,1:n-1) + xSWAN(:,2:n))  xSWAN(:,n)];
% ySWAN = [ ySWAN(:,1)  0.5*(ySWAN(:,1:n-1) + ySWAN(:,2:n))  ySWAN(:,n)];
% xSWAN = [ xSWAN(1,:); 0.5*(xSWAN(1:m-1,:) + xSWAN(2:m,:)); xSWAN(m,:)];
% ySWAN = [ ySWAN(1,:); 0.5*(ySWAN(1:m-1,:) + ySWAN(2:m,:)); ySWAN(m,:)];
% clear m n
%
% [m n] = size(lonC);
% xC = lonC; yC = latC;
% xC = [ xC(:,1)  0.5*(xC(:,1:n-1) + xC(:,2:n))  xC(:,n)];
% yC = [ yC(:,1)  0.5*(yC(:,1:n-1) + yC(:,2:n))  yC(:,n)];
% xC = [ xC(1,:); 0.5*(xC(1:m-1,:) + xC(2:m,:)); xC(m,:)];
% yC = [ yC(1,:); 0.5*(yC(1:m-1,:) + yC(2:m,:)); yC(m,:)];
% clear m n
% rmpath('C:\Users\sdalyander\Documents\MATLAB\m_cmg\trunk\tri\')
% SWANboxes = NaN(size(locs)); Cboxes = NaN(size(locs));
% isGood = find(maskSWAN(:) == 1);
% for bb = 1:length(isGood)
%     [ii,jj] = ind2sub(size(lonSWAN),isGood(bb));
%     nodes = [xSWAN(ii,jj) ySWAN(ii,jj); xSWAN(ii,jj+1) ySWAN(ii,jj+1); ...
%         xSWAN(ii+1,jj+1) ySWAN(ii+1,jj+1); xSWAN(ii+1,jj) ySWAN(ii+1,jj)];
%     [isIn,isOn] = inpoly(locs, nodes);
%     isIn = isIn | isOn; clear isOn
%     SWANboxes(isIn,1) = ii;
%     SWANboxes(isIn,2) = jj;
% end
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
nc = mDataset('dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CompGrid/swan_agg.ncml');

zrMax = 5;
ncC = mDataset(mp.urlC);
fTime = value2Index(timeC,min(timeSWAN));
lTime = value2Index(timeC,max(timeSWAN));
done = zeros(mp.NL,1);

Cboxes(Cboxes(:,1) == size(h,1),:) = NaN;
Cboxes(Cboxes(:,2) == size(h,2),:) = NaN;

%% 
for nn2 = 987:mp.NL
    disp(['On n = ' num2str(nn2) ' of ' num2str(mp.NL) '.'])
    if done(nn2) || isnan(SWANboxes(nn2,1)) || isnan(Cboxes(nn2,1)), continue, end
    
    zetaA = double(ncC{'zeta'}(fTime:lTime,Cboxes(nn2,1),Cboxes(nn2,2)));
    
    %Pull the wave data
    disp('Pulling down data')
    tic
    ubrA = nc{'ub'}(:,SWANboxes(nn2,1),SWANboxes(nn2,2));
    tbrA = nc{'tmbot'}(:,SWANboxes(nn2,1),SWANboxes(nn2,2));
    wdirA = nc{'bdir'}(:,SWANboxes(nn2,1),SWANboxes(nn2,2));
    wdirA = (wdirA*pi/180)*-1 + (3*pi/2); %Degrees relative to North, direction waves coming from to radians rel. East, direction to
    
    ll = 1;
    u = double(squeeze(ncC{'u'}(fTime:lTime,ll,Cboxes(nn2,1),Cboxes(nn2,2)-1:Cboxes(nn2,2))));
    v = double(squeeze(ncC{'v'}(fTime:lTime,ll,Cboxes(nn2,1)-1:Cboxes(nn2,1),Cboxes(nn2,2))));
    currA = nanmean(u,2) + 1i.*nanmean(v,2); clear u v
    toc
    
    isThis = find(SWANinds == SWANinds(nn2) & Cinds == Cinds(nn2));
    
    for nn3 = 1:length(isThis)
        nn = isThis(nn3);
        zrN = NaN(length(timeSWAN),1);
        tauwc = NaN(length(timeSWAN),1);
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
        
        disp('Starting time analysis')
        tic
        for tt = 1:length(timeSWAN)
            %Check current model depth
            
            zeta = interp1(timeC(fTime:lTime),zetaA,timeSWAN(tt));
            z_r = z_rB(:,sub2ind(size(h),Cboxes(nn,1),Cboxes(nn,2))) + ...
                scmCshc*zeta/h(Cboxes(nn,1),Cboxes(nn,2))' + ...
                (1+Cs_r)*zeta';
            z_w = z_wB(:,sub2ind(size(h),Cboxes(nn,1),Cboxes(nn,2))) + ...
                scmCshc_w*zeta/h(Cboxes(nn,1),Cboxes(nn,2))' + ...
                (1+Cs_w)*zeta';
            zrA = z_r - repmat(z_w(1),[length(z_r) 1]);
            
            if zrA(1) > zrMax, continue, end
            
            
            curr = interp1(timeC(fTime:lTime),currA,timeSWAN(tt));
            
            Umag = sqrt(real(curr).^2 + imag(curr).^2 + eps);
            phic = atan2(imag(curr),real(curr)) + angle_r(Cboxes(nn,1),Cboxes(nn,2));
            phiwc = abs(wdirA(tt) - phic);
            zr = zrA(ll);
            
            m = m94(ubrA(tt),2*pi/tbrA(tt),Umag,zr,phiwc,kN,0);
            
            
            while m.dwc > zr
                ll = ll+1;
                gg = value2Index(timeC, timeSWAN(tt));
                u2 = double(squeeze(ncC{'u'}(gg-1:gg+1,ll,Cboxes(nn,1),Cboxes(nn,2)-1:Cboxes(nn,2))));
                v2 = double(squeeze(ncC{'v'}(gg-1:gg+1,ll,Cboxes(nn,1)-1:Cboxes(nn,1),Cboxes(nn,2))));
                curr = nanmean(u2,2) + 1i.*nanmean(v2,2);
                curr = interp1(timeC(gg-1:gg+1),curr,timeSWAN(tt));
                
                clear u v
                Umag = sqrt(real(curr).^2 + imag(curr).^2 + eps);
                phic = atan2(imag(curr),real(curr)) + angle_r(Cboxes(nn,1),Cboxes(nn,2));
                phiwc = abs(wdirA(tt) - phic);
                zr = zrA(ll);
                
                m = m94(ubrA(tt),2*pi/tbrA(tt),Umag,zr,phiwc,kN,0);
            end %Ends while loop checking zr
            zrN(tt) = zr;
            tauwc(tt) = 1025*m.ustrr.^2;
        end %Ends time series loop
        toc
        
        for ii = 1:length(mmList)
            
            myTime =  (gTime(:,1) == yyList(ii)) & (gTime(:,2) == mmList(ii));
            
            ncR = netcdf(fileN{ii},'write');
            ncR{'lon'}(nn) = locs(nn,1);
            ncR{'lat'}(nn) = locs(nn,2);
            ncR{'kN'}(nn) = kN;
            ncR{'tau_crit'}(nn) = tau_crit;
            ncR{'zr'}(:,nn) = zrN(myTime);
            ncR{'tauwc'}(:,nn) = tauwc(myTime);
            close(ncR)
        end
    end %Ends loop of points with same boxes
    
end
close(nc), clear nc
close(ncC), clear ncC
