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
load analysis3_spec_welch
mVar = squeeze(spec_analysis_low(1,1,:,:)./spec_analysis_high(1,1,:,:));
freqStressRaw=mVar;
% figure, subplot(1,2,1)
% pcolorjw(lon,lat,freqStressRaw)
% title('Frequency of Stress')
cuts = [-Inf 2 10 Inf];
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

%% Magnitude of Disturbance
% % What % of the bed is moved by median (or 90th percentile) stress
% load perc_grainsize_crit_smaller_percentiles
load analysis1_gen_geo lon lat mmean
test = squeeze(mmean(1,1,:,:));
clear mmean
% test(~isnan(test)) = 1;
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
load compare_sediment_percentiles_to_skinstress_tseries_cohesive
thresh = 20;
mc = obs_model_percentiles(:,1);
p = locs;

goods = find(~isnan(mc));
F = TriScatteredInterp(double(p(goods,1)),double(p(goods,2)),double(mc(goods)));
MC = F(double(lon),double(lat)).*test;

% figure, subplot(1,2,1), pcolorjw(lon,lat,MC)
% title('Frequency of Disturbance')
mVar = MC;
freqDisturbRaw = mVar;
cuts = [-Inf 5 50 Inf];
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

% 
% %% manual analysis
% % close all
% figure
% mSort = [magStress(:) freqStress(:) magDisturb(:) freqDisturb(:)];
% % mSort = cidx2;
% [outMat,count] = groupbypermute(mSort(:,1:4));
% 
% %Get rid of small ones
% outMatNew = outMat;
% % bads = find(count(:,2)<=5);
% % outMatNew(ismember(outMatNew,count(bads,1))) = NaN;
% subplot(1,2,1),pcolorjw(lon,lat,reshape(outMatNew,size(lon)))
% title({[num2str(size(count,1)) ' Regions out of']; ...
%     [num2str(3^size(mSort,2)) ' Possible Combinations']})
% dasp(38.5)
% %%
% [sortList,ii] = sort(count(:,2),'descend');
% cumList = 100*cumsum(sortList)/length(find(~isnan(outMat)));
% cut = find(cumList > 95,1,'first')
% 
% %%
% 
% % figure
% % jj = jj+1
% % mSort(find(outMat(:) == count(ii(jj),1),1,'first'),1:4)
% % outMatNew = outMat; outMatNew(~isnan(outMatNew)) = 0; outMatNew(outMat == count(ii(jj),1)) = 1;
% % figure, subplot(1,2,1),pcolorjw(lon,lat,reshape(outMatNew,size(lon)))
% % dasp(38.5)
% %1. Southern MAB, mid-shelf (2,3,3,2)*
% %2. Offshore & Western MP (1,2,1,1)**
% %3. Northern MAB, mid-shelf (2,2,3,2)*
% %4. Offshore (1,1,1,1) &&
% %5.  Northern MAB, mid-shelf band (1,2,2,1)**
% %6. Eastern MP, south HSV (1,2,3,2)^
% %7.  Inshore near Delaware Bay (3,2,3,2)***
% %8.  Mid-Mud Patch (1,2,2,2)@
% %9.  Southern MAB, nearshore (3,3,3,2)***
% %10.  Northern MAB, mid-shelf (1,2,3,1)^
% %11. Mid MAB, mid-shelf (2,2,3,1)^
% %12.  Nantucket Shoals (3,1,3,3)^^
% %13.  south Nantucket Shoals (2,1,3,2)^^
% %14.  In the bays (1,1,3,3) &
% %15.  Bay mouths, edge Nantucket Shoals (3,2,3,3)^^
% %16.  Offshore of Nantucket Shoals (1,1,3,2) &
% %17.  Offshore and bay pieces (1,1,2,2) &&
% %18.  mid-shelf pieces (1,3,2,1)**
% %19.  Bays (2,1,3,3) &
% %20.  Mid-shelf pieces (1,3,3,1)%%
% %21.  Mid-shelf pieices (2,3,3,1)%%
% %22.  NE Nantucket Shoals (3,1,2,3)^^
% %23.  southern MAB mid shelf pieces (1,3,3,2)*
% %24.  norther MAB offshore pieces (1,1,2,1) &&
% %25.  band around Nantucket Shoals (3,1,3,2)^^
% %26. mid-shelf blobs (2,2,2,1)**
% %27.  Over Nantucket Shoals (3,1,2,2)^^
% %28.  scattered pieces (1,2,3,3)
% %29.  mid-shelf blob (2,3,2,2)^^
% %30.  nearshore pieces (3,3,3,3)***
% %31. mid-shelf pieces (2,2,2,2)@
% %32.  Bays pieces (2,1,2,2)^^
% %33.  Bits (2,2,3,3)
% %34.  Offshore bits (1,1,1,2)&&
% %35.  Offshore bits (1,1,3,1)%%
% %36.  Midshelf bits (1,3,1,1)**
% %37.  Offshore bits (1,3,2,2)
% %38.  Bits (2,3,2,1) **
% %39. Bits (1,2,1,2)**
% %40. Bits (2,3,3,3)
% %41.  Bits (1,1,2,3) ^^
% %42.  Bits (2,1,1,1)
% %43.  Bits (3,1,1,1)
% %44.  Bits (3,1,1,2)
% %45.  Bits (3,2,2,2)^^
% %46.  Bits (2,1,2,1)
% %47.  Bits (3,3,2,2)^^
% %48. Bits (3,1,1,3)
% %49.  Bits (3,2,2,1)
% 
% %%
% outMatNew = NaN(size(outMat));
% 
% %Region 1:  Low to moderate Stress, storm or mixed dominated, deep
% %disturbance occuring infrequently or moderately often
% outMatNew(ismember(outMat,[count(ii(1),1) count(ii(3),1) ...
%     count(ii(23),1) count(ii(6),1) ...
%     count(ii(10),1) count(ii(11),1)]))=1;
% %Region 2: Low to moderate stress, storm or mixed dominated, low to moderate disturbance occuring
% %infrequently
% outMatNew(ismember(outMat,[count(ii(2),1) count(ii(5),1) ....
%     count(ii(18),1) count(ii(26),1) count(ii(36),1) count(ii(38),1)...
%     count(ii(39),1)])) = 2;
% %Region 3:  Low stress, tidal dominated, low or moderate disturbance occurring
% %infrequently or moderately
% outMatNew(ismember(outMat,[count(ii(4),1) count(ii(24,1)) ...
%     count(ii(17),1) count(ii(34),1)])) = 3;
% %Region 4: High stress, mixed or storm dominated, deep disturbance occuring
% %moderately or highly often
% outMatNew(ismember(outMat,[count(ii(7),1) count(ii(9),1) count(ii(30),1)])) = 4;
% %Region 5:High stress, mixed or tidally dominated, moderate to high disturbance occurring
% %moderately to frequently
% outMatNew(ismember(outMat,[count(ii(12),1) count(ii(15),1) ...
%     count(ii(22),1) count(ii(25),1) count(ii(27),1) count(ii(45),1) count(ii(47),1)])) = 5;
% %Region 6:Low to Moderate stress, tidally dominated, moderate to high disturbance occurring
% %moderately to frequently
% outMatNew(ismember(outMat,[count(ii(13),1) count(ii(14),1) count(ii(19),1) ...
%     count(ii(16),1) count(ii(29),1) count(ii(32),1) count(ii(41),1)])) = 6;
% %Region 7:  Low to moderate Stress, mixed dominated, moderate disturbance occurring
% %moderately frequently
% outMatNew(ismember(outMat,[count(ii(8),1) count(ii(31),1)])) = 7;
% %Region 8:  Low to moderate stress, tidally dominated, high disturbance
% %occurring infrequently
% outMatNew(ismember(outMat,[count(ii(20),1) count(ii(21),1) count(ii(35),1)])) = 8;
% %Region 9:  Not in a zone
% % outMatNew(find(~isnan(outMat) & isnan(outMatNew))) = 9;
% 
% subplot(1,2,2), pcolorjw(lon,lat,reshape(outMatNew,size(lon))), colorbar
% title('Distilled to 8 Groups (Some Points Not Assigned)')
% dasp(38.5)
