%esp_plots

xL = [-81.5 -74];
yL = [24 37.5];

xfac=cos(28*pi/180);
lab = get(cbar,'xlabel');
set(cbar,'fontname','times','fontsize',14)
set(lab,'string',cBarLabel,'fontsize',14,'fontname','times')

thePos2 = get(mAx,'position');
thePos = thePos2;
thePos(4) = thePos(4)*1.2;
thePos(2) = thePos(2) - thePos2(4)*1.2/4+0.09;
set(mAx,'position',thePos)
hold on, box on
% xL = xlim;
% yL = ylim;

ax2 = axes('position',thePos);
hold on
set(ax2,'xlim',xL,'ylim',yL)
cList = -1*[20:20:200];
% cList = [20:20:200];
[c2,h] = contour(lon1,lat1,topo,[cList(:) cList(:)]);
% [c2,h] = contour(lonAll-360,latAll,-1*depthAll,[cList(:) cList(:)]);
% [c2v2,hv2] = contour(lon1,lat1,topo,[cList2(:) cList2(:)]);
set(h,'color','k','linewidth',0.2);
% set(hv2,'color',[0.3 0.3 0.3],'linewidth',0.2)
set(ax2,'color','none','xtick',[],'ytick',[])
leg = legend(h,'20 m contours');
set(leg,'fontsize',14,'box','off','fontname','times','position',[0.39 0.5 0.2964 0.0675])
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

%Scale bar...at 28.5 degrees latitude, one degree of longitude = 68.83 miles
%one degree longitude = 97.90 km
mileFac = 68.83;
kmFac = 97.90;
sPoint = -78.4;
yPoint = 28.5;
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
set(mAx,'xlim',xL)
set(mAx,'ylim',yL)
set(mAx,'ytick',[],'xtick',[],'box','off')
if exist('mAx2','var')
set(mAx2,'position',thePos)
set(mAx2,'xlim',xL)
set(mAx2,'ylim',yL)
set(mAx2,'ytick',[],'xtick',[],'box','off')
end
g = get(cbar,'position');
% set(cbar,'dataaspect',[25 1 1]);
g2 = g;
g2(3) = g(3)/4;
g2(1) = (g(1)+g(3)/2) - g2(3)/2;
 g2(2) = g(2) - 0.04;
g2(4) = g2(4)/2;
set(cbar,'position',g2)
