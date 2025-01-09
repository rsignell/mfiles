%esp_plots

xfac=cos(38*pi/180);
lab = get(cbar,'xlabel');
set(cbar,'fontname','times','fontsize',14)
set(lab,'string',cBarLabel,'fontsize',14,'fontname','times')

thePos2 = get(mAx,'position');
thePos = thePos2;
thePos(4) = thePos(4)*1.2;
thePos(2) = thePos(2) - thePos2(4)*1.2/4+0.09;
set(mAx,'position',thePos)
hold on, box on
xL = xlim;
yL = ylim;
ax2 = axes('position',thePos);
hold on
set(ax2,'xlim',xL,'ylim',yL)
cList = -1*[20:20:100 200:50:400];
[c2,h] = contour(lon1,lat1,topo,[cList(:) cList(:)]);
set(h,'color','k','linewidth',0.2);
set(ax2,'color','none','xtick',[],'ytick',[])
leg = legend(h,'20 m contours','location','southeast');
set(leg,'fontsize',14,'box','off','fontname','times')
[c3,h2] = contourf(lon1,lat1,topo,[0 0]);
ch = get(h2,'children');
for ii = 1:length(ch)
set(ch(ii),'facecolor',[0.9 0.9 0.9])
end
axis([xL yL])
for ii = 1:length(PB.seg)
    loc = PB.pts{ii};
    if any((loc(:,1) > min(ylim)) & (loc(:,1) < max(ylim)) & ...
            (loc(:,2) > min(xlim)) & loc(:,2) < max(xlim));
        plot(loc(:,2),loc(:,1),'k')
    end
end
for ii = 1:length(PB2.seg)
    loc = PB2.pts{ii};
    if any((loc(:,1) > min(ylim)) & (loc(:,1) < max(ylim)) & ...
            (loc(:,2) > min(xlim)) & loc(:,2) < max(xlim));
        plot(loc(:,2),loc(:,1),'k')
    end
end
%Added 5/10/2011
xlim(xL);
ylim(yL);

%Scale bar...at 35 degrees latitude, one degree of longitude = 56.72 miles
%one degree longitude = 91.29 km
mileFac = 56.72;
kmFac = 91.29;
sPoint = -70.7;
yPoint = 35;
nMiles = 100;
nKm = 100;
l = line([sPoint sPoint+nMiles/mileFac],[yPoint yPoint]);
set(l,'linewidth',1.5,'color','k');
t = text(min([sPoint sPoint+nMiles/mileFac]),yPoint,[num2str(nMiles) ' miles']);
set(t,'fontname','times','fontsize',14,'horiz','left','vert','bottom')

nudge = 0.1;
l = line([sPoint sPoint+nKm/kmFac],[yPoint-nudge yPoint-nudge]);
set(l,'linewidth',1.5,'color','k');
t = text(min([sPoint sPoint+nKm/kmFac]),yPoint-nudge,[num2str(nMiles) ' km']);
set(t,'fontname','times','fontsize',14,'horiz','left','vert','top')

set(gca,'xticklabel','','yticklabel','','DataAspectRatio',[1 xfac 1000])
set(mAx,'xticklabel','','yticklabel','','DataAspectRatio',[1 xfac 1000])
if exist('mAx2','var')
set(mAx2,'xticklabel','','yticklabel','','DataAspectRatio',[1 xfac 1000])
end
set(ax2,'position',thePos)
set(mAx,'position',thePos)
if exist('mAx2','var')
set(mAx2,'position',thePos)
end
g = get(cbar,'position');
% set(cbar,'dataaspect',[25 1 1]);
g2 = g;
g2(3) = g(3)/4;
g2(1) = (g(1)+g(3)/2) - g2(3)/2;
 g2(2) = g(2) - 0.04;
g2(4) = g2(4)/2;
set(cbar,'position',g2)
