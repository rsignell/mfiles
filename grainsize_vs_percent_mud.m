%%grainsize_vs_percent_mud

ncclear, clc
load perc_grainsize_crit_smaller_percentiles mc p gsizel

split = find(gsizel < 0.004e-3,1,'last');

percMud = NaN(size(mc,1),1);
d90 = NaN(size(mc,1),1);
d50 = NaN(size(mc,1),1);

for gg = 1:size(mc,1)
    disp(['On gg = ' num2str(gg) ' of ' num2str(size(mc,1)) '.'])
    gall = mc(gg,:); gall(isnan(gall)) = 0;
    if all(gall==0), continue, end
    gind = gall;
    gall = cumsum(gall);
    
    tt = find(gall >= 90,1,'first');
    d90(gg) = gsizel(tt);
    tt = find(gall >= 50,1,'first');
    d50(gg) = gsizel(tt);
    
    percMud(gg) = gall(split);
end

isCohesive = find(percMud >= 7.5);