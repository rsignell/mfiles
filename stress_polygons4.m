%stress_polygons

%Want to set up indices to characterize:
%1.  Magnitude of stress (mean or 95%)
%2.  Frequency of stress (energy in high vs. low frequency)
%3.  Disturbance

ncclear, close all

%% Magnitude of stress
close all
load analysis1_gen_geo m95 lon lat
mVar = squeeze(m95(1,1,:,:));
magStressRaw = mVar;
cuts = [-Inf 0.5 1 Inf];
% figure, subplot(1,2,1)
% pcolorjw(lon,lat,magStressRaw)
% title('Magntitude of Stress')
polysort = NaN(size(mVar));
for cc = 1:length(cuts)-1
    polysort(mVar >= cuts(cc) & mVar < cuts(cc+1)) = cc;
end
clear cc
subplot(2,2,1)
pcolorjw(lon,lat,polysort)
dasp(38.5)
magStress = polysort;
title({'Magnitude of Stress (95th Percentile)'; ['< ' num2str(cuts(2)) ', ' num2str(cuts(2)) '-' num2str(cuts(3)) ', ' ...
    '> ' num2str(cuts(3))]})

%% Frequency of stress
load analysis3_spec_welch spec_analysis
mVar = squeeze(spec_analysis(1,1,:,:));
freqStressRaw=mVar;
% figure, subplot(1,2,1)
% pcolorjw(lon,lat,freqStressRaw)
% title('Frequency of Stress')
cuts = [-Inf 1 2 Inf];
polysort = NaN(size(mVar));
for cc = 1:length(cuts)-1
    polysort(mVar >= cuts(cc) & mVar < cuts(cc+1)) = cc;
end
clear cc
subplot(2,2,2)
pcolorjw(lon,lat,polysort)
dasp(38.5)
freqStress = polysort;
title({'Low Frequency/High Frequency'; ['< ' num2str(cuts(2)) ', ' num2str(cuts(2)) '-' num2str(cuts(3)) ', ' ...
    '> ' num2str(cuts(3))]})

% %% Magnitude of Disturbance
% % What % of the bed is moved by median (or 90th percentile) stress
% load perc_grainsize_crit_smaller_percentiles
load analysis1_gen_geo lon lat mmean
test = squeeze(mmean(1,1,:,:));
clear mmean
test(~isnan(test)) = 1;
% thresh = 90;
% mc = perc_crit_smaller_obs_all(:,thresh/10);
% goods = find(~isnan(mc));
% F = TriScatteredInterp(p(goods,1),p(goods,2),mc(goods));
% MC = F(lon,lat).*test;
% %  figure, subplot(1,2,1), pcolorjw(lon,lat,MC)
% % title('Magnitude of Disturbance')
% mVar = MC;
% magDisturbRaw = mVar;
% cuts = [-Inf 25 75 Inf];
% polysort = NaN(size(mVar));
% for cc = 1:length(cuts)-1
%     polysort(mVar >= cuts(cc) & mVar < cuts(cc+1)) = cc;
% end
% clear cc
% subplot(2,2,3)
% pcolorjw(lon,lat,polysort)
% dasp(38.5)
% magDisturb = polysort;
% title({['% Bed Mobilized by ' num2str(thresh) 'th Stress']; ['< ' num2str(cuts(2)) ', ' num2str(cuts(2)) '-' num2str(cuts(3)) ', ' ...
%     '> ' num2str(cuts(3))]})


%% Frequency of Disturbance
% How often is 20%+ of the bed moved?
load ../point_stress/compare_sediment_percentiles_to_skinstress_tseries_cohesive
thresh = 20;
mc = obs_model_percentiles(:,1);
p = locs;
goods = find(~isnan(mc));
F = TriScatteredInterp(p(goods,1),p(goods,2),mc(goods));
MC = F(lon,lat).*test;

% figure, subplot(1,2,1), pcolorjw(lon,lat,MC)
% title('Frequency of Disturbance')
mVar = MC;
freqDisturbRaw = mVar;
cuts = [-Inf 25 75 Inf];
polysort = NaN(size(mVar));
for cc = 1:length(cuts)-1
    polysort(mVar >= cuts(cc) & mVar < cuts(cc+1)) = cc;
end
clear cc
subplot(2,2,4)
pcolorjw(lon,lat,polysort)
dasp(38.5)
freqDisturb = polysort;
title({['% of Time ' num2str(thresh) '% of Bed is Mobilized']; ['< ' num2str(cuts(2)) ', ' num2str(cuts(2)) '-' num2str(cuts(3)) ', ' ...
    '> ' num2str(cuts(3))]})

%% Correlation waves/currents
% load analysis3_corr_spec xcor_wave_curr
% mVar = squeeze(xcor_wave_curr(1,:,:));
% corrStressRaw=mVar;
% figure, subplot(1,2,1)
% pcolorjw(lon,lat,corrStressRaw)
% title('Correlation Waves/Current')
% cuts = [-Inf 0.1 0.5 Inf];
% polysort = NaN(size(mVar));
% for cc = 1:length(cuts)-1
%     polysort(mVar >= cuts(cc) & mVar < cuts(cc+1)) = cc;
% end
% clear cc
% subplot(1,2,2)
% pcolorjw(lon,lat,polysort)
% corrStress = polysort;
% title({'Correlation Waves/Currents'; ['< ' num2str(cuts(2)) ', ' num2str(cuts(2)) '-' num2str(cuts(3)) ', ' ...
%     '> ' num2str(cuts(3))]})



%% Kmeans
% % close all
% nVars = 3;
% figure
% clear cidx2 cmeans2
% mVars = [magStressRaw(:) freqStressRaw(:) magDisturbRaw(:) freqDisturbRaw(:) corrStressRaw(:)];
% titles = {'Magnitude of Stress'; 'Frequency of Stress'; 'Magnitude of Disturbance'; ...
%     'Frequency of Disturbance'; 'Strength of Wind-Driven Currents'};
% for ss = 1:5
%     subplot(2,3,ss)
% [cidx2(:,ss),cmeans2(:,ss)] = kmeans(mVars(:,ss),nVars);
% pcolorjw(lon,lat,reshape(cidx2(:,ss),size(lon)))
% title(titles{ss})
% dasp(41)
% end

%% Groups
% close(8)
% ii = 5
% [idx,cmeans] = kmeans(cidx2(:,[1:5]),ii,'distance','sqEuclidean','start','uniform');
% figure, pcolorjw(lon,lat,reshape(idx,size(lon))), title([num2str(ii) ' Groups'])


%% manual analysis
% close all
figure
% mSort = [magStress(:) freqStress(:) magDisturb(:) freqDisturb(:)];
% mSort = [magStress(:) freqStress(:) magDisturb(:)];
mSort = [magStress(:) magDisturb(:) freqDisturb(:)];
% mSort = cidx2;
[outMat,count] = groupbypermute(mSort);

%Get rid of small ones
outMatNew = outMat;
% bads = find(count(:,2)<=5);
% outMatNew(ismember(outMatNew,count(bads,1))) = NaN;
% subplot(1,2,1),pcolorjw(lon,lat,reshape(outMatNew,size(lon)))
figure
% [c,h] = contourf(lon,lat,reshape(outMatNew,size(lon)));
pcolorjw(lon,lat,reshape(outMatNew,size(lon)))
title({[num2str(size(count,1)) ' Regions out of']; ...
    [num2str(3^size(mSort,2)) ' Possible Combinations']})
text(-73,37,{'Criteria:'; 'Mag. Stress,'; 'Mag. Dist.,'; 'Freq. Dist.'});

dasp(38.5)
%%
[sortList,ii] = sort(count(:,2),'descend');
cumList = 100*cumsum(sortList)/length(find(~isnan(outMat)));
cut = find(cumList > 95,1,'first')

%%

% close
% jj = 10
% mSort(find(outMat(:) == count(ii(jj),1),1,'first'),:)
% outMatNew = outMat; outMatNew(~isnan(outMatNew)) = 0; outMatNew(outMat == count(ii(jj),1)) = 1;
% figure, pcolorjw(lon,lat,reshape(outMatNew,size(lon)))


%1. Mid-shelf band (2,3,2) 3
%2. Off-shore band (1,1,1) 1
%3. Outer mid-shelf band (1,2,1) 6
%4. Nearshore band (3,3,2) 3
%5. West Mud Patch around Nantucket Shoals (1,3,2) 5
%6. West Mud Patch around Nantucket Shoals (1,2,2) 5
%7. Nantucket Shoals (3,3,3) 2 
%8.  Nearshore band (1,3,1) 6
%9.  Mid-shelf band (2,3,1) 4
%10.  Bits (1,3,3) 5
%11.  Bits (2,3,3) 2
%12.  Bits (2,2,2)3
%13.  Bits (2,2,1) 4
%14.  Nantucket Shoals (3,2,3) 2 
%15.  Nantucket Shoals (3,2,2) 2
%16.  Offshore bits (1,1,2) 1
%17.  Bay bits (1,2,3) 5
%18.  Bits (2,1,1) 1
%19.  Nantucket Shoals (3,1,1) 1
%20.  Nantucket Shoals (3,1,2) 1
%21.  Nantucket Shoals (3,1,3) 1
%22.  Bits(3,2,1) 4

outMatNew = NaN(size(outMat));
%Region 1:  Any Stress, Low Disturbance at any frequency
outMatNew(ismember(outMat,count(ii([2 16 18 19 20 21]),1))) = 1;

%Region 2:  Moderate-high Stress, Moderate-high Disturbance, Tides
outMatNew(ismember(outMat,count(ii([7 11 14 15]),1))) = 2;

%Region 3:  Moderate-high Stress, Moderate-high Disturbance, Mixed
outMatNew(ismember(outMat,count(ii([1 4 12]),1))) = 3;

%Region 4:  Moderate-high Stress, Moderate-high Disturbance, Storms
outMatNew(ismember(outMat,count(ii([9 13 22]),1))) = 4;

%Region 5:  Low Stress, Moderate-high disturbance, Mixed-Tides
outMatNew(ismember(outMat,count(ii([5 6 10 17]),1))) = 5;

%Region 6:  Low Stress, Moderate-high Disturbance, Storms
outMatNew(ismember(outMat,count(ii([3 8]),1))) = 6;


% subplot(1,2,2), pcolorjw(lon,lat,reshape(outMatNew,size(lon)))%, colorbar
figure
% [c,h] = contourf(lon,lat,reshape(outMatNew,size(lon)));
pcolorjw(lon,lat,reshape(outMatNew,size(lon)))
% clabel(c,h)
title('Distilled to 6 Groups')
c = colorbar('west');
set(c,'ytick',[1:1:6])
dasp(38.5)
