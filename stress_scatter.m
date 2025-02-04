%% Stress scatter

load analysis1_gen_geo.mat      
load analysis4_disturbance.mat

mp.urlC = ...
   'http://tashtego.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2009_da/his';
nc = mDataset(mp.urlC);
h = nc{'h'}(:);
close(nc); clear nc



%% Plot
close all
hs = find(h(:) > 200);
testx = squeeze(mmed(1,1,:,:));
testy = crit_stress;
testy2 = mean_gsize;
labx = 'Median Stress (m)';
laby = 'Critical Stress (Pa)';
laby2 = 'Grain Size (m)';
limsy = [0 2];
limsy2 = [0 0.001];
limsx = [0 2];
dps = [0:max(limsx)/100:max(limsx)];
figure, hold on
subplot(1,2,1), hold on
xVar = testx(:);
yVar = testy(:);
plot(xVar,yVar,'k.')
% plot(xVar(hs),yVar(hs),'r.')
axis equal
axis([limsx limsy])
xlabel(labx)
ylabel(laby)
l = line(limsx, limsy);
set(l,'color','k','linestyle','-.')
good = find(~isnan(xVar + yVar));
pf = polyfitn(xVar(good),yVar(good),1);
l2 = plot(dps, polyvaln(pf,dps));
res = yVar - polyvaln(pf,xVar);
res = reshape(res,size(lon));
R = corrcoef(xVar(good),yVar(good));
set(l2,'color','r')
box on

subplot(1,2,2), hold on
xVar = testx(:);
yVar = testy2(:);
plot(xVar,yVar,'k.')
% axis equal
axis([limsx limsy2])
xlabel(labx)
ylabel(laby2)
l = line(limsx, limsy2);
set(l,'color','k','linestyle','-.')
good = find(~isnan(xVar + yVar));
pf2 = polyfitn(xVar(good),yVar(good),1);
l2 = plot(dps, polyvaln(pf,dps));
res2 = yVar - polyvaln(pf,xVar);
res2 = reshape(res,size(lon));
R2 = corrcoef(xVar(good),yVar(good));
set(l2,'color','r')
box on
set(gca,'plotboxaspectratio',[1.1 1 1])