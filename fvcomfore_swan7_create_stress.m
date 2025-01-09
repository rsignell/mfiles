ncclear

url = 'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_CFSR_7grid_v2/swan_agg.ncml';

nc = ncgeodataset(url);
lon = nc.data('lon');
lat = nc.data('lat');
time = nc.data('time') + datenum(1858,11,17,00,00,00);
gTime = datevec(time);

yrs = [2010 2010 2010 2011 2011 2011 2011];
mos = [10 11 12 1 2 3 4];

fLon = min(find(lon >= -74+360));
lLon = length(lon);
fLat = min(find(lat >= 40));
lLat = length(lat);

depth = nc.data('depth',[fLat fLon],[lLat lLon]);

mp.history = 'FVCOM Forecast Currents, SWAN 7 grid CFSR';
mp.type = 'Tau';

%Can either be single value or spatially variant
kN = 0.5e-2;    %5 cm physical roughness
rhoSW = 1025;   %1025 kg/m3 density seawater

lon = lon(fLon:lLon)-360;
lat = lat(fLat:lLat);

[Lon,Lat] = meshgrid(lon,lat);
[mp.RL,mp.CL] = size(Lon);

kN = repmat(kN,size(Lon));

ncF = mDataset('http://fvcom.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_gom3');
timeFVCOM = nj_time(ncF,'u');
lonFVCOM = double(ncF{'lon'}(:));
latFVCOM = double(ncF{'lat'}(:));
loncFVCOM = double(ncF{'lonc'}(:));
latcFVCOM = double(ncF{'latc'}(:));
nv = double(ncF{'nv'}(:));

ind = NaN(1,numel(Lon));
indc = NaN(1,numel(Lon));

% for ii = 1:numel(Lon)
%     ind(ii) = nearxy(lonFVCOM(:),latFVCOM(:),Lon(ii),Lat(ii));
%     indc(ii) = nearxy(loncFVCOM(:),latcFVCOM(:),Lon(ii),Lat(ii));
% end
load('C:\Users\sdalyander\Documents\StressAnalysis\FVCOM_SWAN\FVCOM_inds.mat')
zrMax = 10;
ll_list = NaN(numel(Lon),1);
for tt2 = 2%:length(yrs)
    disp(['Processing tt2 = ' num2str(tt2) ' of ' num2str(length(yrs))])
    
    fname = ['FVCOM_SwanH_kNpt5_' num2str(yrs(tt2)-2000) ...
        num2str(mos(tt2),'%02.0f') '.nc'];
    
    
    fTime = min(find(gTime(:,1) == yrs(tt2) & gTime(:,2) == mos(tt2)));
    lTime = max(find(gTime(:,1) == yrs(tt2) & gTime(:,2) == mos(tt2)));
    fInd = [fTime fLat fLon];
    lInd = [lTime lLat lLon];
    
    disp('Loading wave results')
    ub = nc.data('ub',fInd,lInd);
    tmbot = nc.data('tmbot',fInd,lInd);
    bdir = nc.data('bdir',fInd,lInd);
    bdir = (bdir*pi/180)*-1 + (3*pi/2); %Radians relative to east, direction waves going to
    
    mTime = time(fTime:lTime);
    
    inTime = find(timeFVCOM >= min(mTime) & timeFVCOM <= max(mTime));
    timeFVCOMi = timeFVCOM(inTime);
    mp.tL = length(timeFVCOMi);
    
    status = create_tau_output(fname,mp);
    if status == 0
        error(['Failed to create file ' fname])
    end
    
    ncW = netcdf(fname, 'write');
    ncW{'lon'}(:) = Lon;
    ncW{'lat'}(:) = Lat;
    ncW{'kN'}(:) = kN;
    ncW{'rho'}(:) = rhoSW;
    ncW{'time'}(:) = timeFVCOMi - datenum(1858,11,17,00,00,00);
    
    disp('Starting spatial analysis')
    for ii = 1:numel(Lon)
        if any(isnan(ub(:,ii)))
            [ii2,jj2] = ind2sub(size(Lon),ii);
            ncW{'tauw'}(:,ii2,jj2) = NaN(1,length(timeFVCOMi));
            ncW{'tauc'}(:,ii2,jj2) = NaN(1,length(timeFVCOMi));
            ncW{'tauwc'}(:,ii2,jj2) = NaN(1,length(timeFVCOMi));
            continue
        end
        
        disp(['On ii = ' num2str(ii) ' of ' num2str(numel(Lon))])
        gridCorners = nv(:,indc(ii));
        h = double([ncF{'h'}(gridCorners(1)) ncF{'h'}(gridCorners(2)) ncF{'h'}(gridCorners(3))]);
        F = TriScatteredInterp(lonFVCOM(gridCorners),...
            latFVCOM(gridCorners),h(:));
        
        sigvar = double([ncF{'siglay'}(:,gridCorners(1)) ncF{'siglay'}(:,gridCorners(2)) ncF{'siglay'}(:,gridCorners(3))]);
        zeta = double([ncF{'zeta'}(inTime,gridCorners(1)) ncF{'zeta'}(inTime,gridCorners(2)) ncF{'zeta'}(inTime,gridCorners(3))]);
        
        
        %Interpolate wave model results to the circulation model times
        ubI = interp1(mTime,ub(:,ii),timeFVCOMi);
        tbrI = interp1(mTime,tmbot(:,ii),timeFVCOMi);
        wdirI = interp1(mTime,bdir(:,ii),timeFVCOMi,'nearest');
        
        
        zrRef = NaN(length(timeFVCOMi),1);
        wave_current_stress = NaN(length(timeFVCOMi),1);
        wave_stress = NaN(length(timeFVCOMi),1);
        current_stress = NaN(length(timeFVCOMi),1);
        
        
        u = ncF{'u'}(inTime,end,indc(ii));
        v = ncF{'v'}(inTime,end,indc(ii));
        Umag = sqrt(u.^2 + v.^2 + eps);
        phic = atan2(v,u);  %Current angle
        
        
        for tt = 1:length(timeFVCOMi)
            if rem(tt,100) == 0,disp(num2str(tt)),end
            ll = length(sigvar);
            
            phiwc = abs(wdirI(tt) - phic(tt));
            
            zetaNow = zeta(tt,:);
            zr = sigvar.*(repmat(h(:)'+zetaNow,length(sigvar),1));
            zr_bot = repmat(h(:)'+zetaNow,length(sigvar),1)+zr;
            
            %Interpolate to the grid center
            zrThis = zr_bot(end,:);
            F = TriScatteredInterp(lonFVCOM(gridCorners),latFVCOM(gridCorners),zrThis');
            zrCent = F(loncFVCOM(indc(ii)),latcFVCOM(indc(ii)));
            if zrCent > zrMax, continue, end %Over max reference height
            
            %Calculate stress
            if isnan(ubI(tt) || isnan(Umag(tt))), continue, end
            m = m94(ubI(tt),2*pi/tbrI(tt),Umag(tt),zrCent,phiwc,kN(ii),0);
            clear phiwc
            
            while m.dwc > zr_bot(ll)
                ll = ll-1;
                u2 = ncF{'u'}(inTime(tt),ll,indc(ii));
                v2 = ncF{'u'}(inTime(tt),ll,indc(ii));
                Umag2 = sqrt(u2.^2 + v2.^2 + eps);
                phic2 = atan2(v2,u2);  %Current angle
                phiwc2 = abs(wdirI(tt) - phic2); clear phic2 u2 v2
                
                %Interpolate to the grid center
                zrThis = zr_bot(ll,:);
                F = TriScatteredInterp(lonFVCOM(gridCorners),...
                    latFVCOM(gridCorners),zrThis');
                zrCent = F(loncFVCOM(indc(ii)),latcFVCOM(indc(ii)));
                
                m = m94(ubI(tt),2*pi/tbrI(tt),Umag2,zrCent,phiwc2,kN(ii),0);
                clear Umag2 phiwc2
            end %Ends while loop checking zr
            
            zrRef(tt) = zrCent;
            wave_current_stress(tt) =rhoSW*m.ustrr.^2;
            wave_stress(tt) =rhoSW*m.ustrwm.^2;
            current_stress(tt) =rhoSW*m.ustrc.^2;
        end
        ll_list(ii) = ll;
        wave_stress(wave_stress > 1e3) = NaN;
        current_stress(current_stress > 1e3) = NaN;
        wave_current_stress(wave_current_stress > 1e3) = NaN;
        [ii2,jj2] = ind2sub(size(Lon),ii);
        ncW{'tauw'}(:,ii2,jj2) = wave_stress;
        ncW{'tauc'}(:,ii2,jj2) = current_stress;
        ncW{'tauwc'}(:,ii2,jj2) = wave_current_stress;
        clear wave_stress current_stress wave_current_stress u v
        clear zeta ubI tbrI Umag Umag2 zrCent phiwc phiw wDirI zrRef
    end
    
    close(ncW)
end
close(ncF), clear ncF
clear nc


