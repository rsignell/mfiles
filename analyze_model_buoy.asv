%analyze_model_buoy
ncclear
tic
load(fullfile('C:\Users\sdalyander\Documents\Waves\Hindcast', 'buoymodeldata2010.mat'))
mp.mdir = 'C:\Users\sdalyander\Documents\Waves\Hindcast';
gtime = datevec(timeAll);

%Get rid of Texas Tower buoy...adrift 01/2011
ttbuoy = strcmp('44066',buoyList);
jan11 = (gtime(:,1) == 2011 & gtime(:,2)==01);
tauAll(jan11,ttbuoy,1) = NaN;
clear ttbuoy jan11

symList = {'.';'s';'^';'h';'p';'d';'*';'+';'o';'x'};


%stress types, analysis time periods
mlists = {1:12; [12 1 2]; 3:5; 6:8; 9:11};

%Set 1:  basic stats
set1 = {'mmean'; 'mstd'; 'm95'};
com1s = {'nanmean('; 'nanstd('; 'prctile('};
com1e = {');'; ');'; ',95);'};

%Set 2:  geologic thresholds
gsize = [62e-6 2e-3];


%% Set up empties
%Set 1...
mmean = NaN([length(mlists) length(buoyList) 2]);
mstd = NaN([length(mlists) length(buoyList) 2]);
m95 = NaN([length(mlists) length(buoyList) 2]);

%Set 2...[gsize timeframe(x5) spatial]
mbed_mobil = NaN([length(gsize) length(mlists) length(buoyList) 2]);

ttypes = {'B'; 'H'};
tIn = [1 3];

%% Run the loop
for bb = 1:length(buoyList)
    disp(['Data set ' num2str(bb) ' of ' num2str(length(buoyList))])
    
    for tt = 1:2
        
        tau = squeeze(tauAll(:,bb,tIn(tt)));
        
        %Cycle through
        for mm = 1:length(mlists)
            inThis = find(ismember(gtime(:,2),mlists{mm}));
            
            %Set 1:  General stats
            for ss = 1:length(set1)
                eval([set1{ss} '(mm,bb,tt)=' com1s{ss} 'tau' ...
                    '(inThis)' com1e{ss}])
            end
            
            %Set 2:  Geologic
            for gg = 1:length(gsize)
                [mbed_mobil(gg,mm,bb,tt),tau_crit] = ...
                    pmsoulsby(gsize(gg),tau(inThis));
            end
        end
    end
end
% close(nc), clear nc
clear overCheck com1e com1s down_dur event_dur ss tau_crit test tt
clear gg gtime ii iiL inThis it jj jjL jt mm
disp('Be sure to save results!!')
toc

%% Plots..set 1
close all
mName = {'Year'; 'Winter (12,01,02)'; 'Spring (03,04,05)'; 'Summer (06,07,08)'; ...
    'Fall (09,10,11)'};
sName = {'Mean \tau_w'; 'STD \tau_w'; '95% \tau_w'};
dUse = find(ismember(buoyList,{'44070','44020'}));
for ss = 1:length(set1)
    eval(['pLims = [min(' set1{ss} '(:)) max(' set1{ss} '(:))];']) 
    figure, orient landscape
    for mm = 1:length(mlists)
        subplot(2,3,mm), hold on
        eval(['xVar = ' set1{ss} '(mm,:,1);'])
        eval(['yVar = ' set1{ss} '(mm,:,2);'])
        
        for bb = 1:length(buoyList)
            p2 = plot(xVar(bb),yVar(bb),['k' symList{bb}]);
            if ismember(bb,dUse)
            set(p2,'markerface','r')
            else
                set(p2,'markerface','k','color','k')
            end
        end
        
        bads = find(isnan(xVar) | isnan(yVar));
        xVar(dUse) = NaN;
        yVar(dUse) = NaN;
        xVar(bads) = []; yVar(bads) = [];
        pf = polyfitn(xVar,yVar,1);
        l2 = plot(pLims,polyvaln(pf,pLims),'r');
        t{1}= ['Slope = ' num2str(pf.Coefficients(1),'%.2f')];
        t{3} = ['R^2 = ' num2str(pf.R2,'%.2f')];
        t{2} = ['Intercept = ' num2str(pf.Coefficients(2),'%.2f')];
        title([sName{ss} ',' mName{mm}])
        xlabel('Buoy (Pa)')
        ylabel('Model (Pa)')
        axis equal
            axis([min(pLims) max(pLims) min(pLims) max(pLims)])
            box on, l = line([min(pLims) max(pLims)], [min(pLims) max(pLims)]);
        set(l,'color','k','linestyle','-.')
        t2 = text(max(xlim), min(ylim), t);
        set(t2,'horiz','right','vert','bottom','fontsize',8)
        clear t
    end
    subplot(2,3,6), hold on
    clear p2
    for bb = 1:10
        p2(bb) = plot(1,1,['k' symList{bb}]);
        if ismember(bb,dUse)
        set(p2(bb),'markerface','r','visible','off')
        else
            set(p2(bb),'markerface','k','visible','off','color','k')
        end
    end
    legend(p2,buoyList)
    axis off
    print('-dpng',['comp_buoy_' set1{ss} '.png'])
end

%% Plot, set 2 (mbed_mobil)
close all
mName = {'Year'; 'Winter (12,01,02)'; 'Spring (03,04,05)'; 'Summer (06,07,08)'; ...
    'Fall (09,10,11)'};
sName = 'Bed Mobility';
gName = {'62 \mum'; '2 mm'};
gFile = {'62micron'; '2mm'};
dUse = find(~ismember(buoyList,{'44070','44020'}));
for gg = 1:length(gsize)
    test = mbed_mobil(gg,:,:,:);
    pLims = [min(test(:)) max(test(:))];
    figure, orient landscape
    for mm = 1:length(mlists)
        subplot(2,3,mm), hold on
        xVar = squeeze(mbed_mobil(gg,mm,:,1));
        yVar = squeeze(mbed_mobil(gg,mm,:,2));
        
        for bb = 1:length(buoyList)
            p2 = plot(xVar(bb),yVar(bb),['k' symList{bb}]);
            if ismember(bb,dUse)
            set(p2,'markerface','r')
            else
                set(p2,'markerface','k','color','k')
            end
        end
        
        xVar(dUse) = NaN; yVar(dUse) = NaN;
        bads = find(isnan(xVar) | isnan(yVar));
        pf = polyfitn(xVar,yVar,1);
        l2 = plot(pLims,polyvaln(pf,pLims),'r');
        t{1}= ['Slope = ' num2str(pf.Coefficients(1),'%.2f')];
        t{3} = ['R^2 = ' num2str(pf.R2,'%.2f')];
        t{2} = ['Intercept = ' num2str(pf.Coefficients(2),'%.2f')];
        title([sName ', ' gName{gg} ', ' mName{mm}])
        xlabel('Buoy (Pa)')
        ylabel('Model (Pa)')
        axis equal
            axis([min(pLims) max(pLims) min(pLims) max(pLims)])
            box on, l = line([min(pLims) max(pLims)], [min(pLims) max(pLims)]);
        set(l,'color','k','linestyle','-.')
        t2 = text(max(xlim), min(ylim), t);
        set(t2,'horiz','right','vert','bottom','fontsize',8)
        clear t
    end
    subplot(2,3,6), hold on
    clear p2
    for bb = 1:10
        p2(bb) = plot(1,1,['k' symList{bb}]);
        if ismember(bb,dUse)
        set(p2(bb),'markerface','r','visible','off')
        else
            set(p2(bb),'markerface','k','visible','off','color','k')
        end
    end
    legend(p2,buoyList)
    axis off
    print('-dpng',['comp_buoy_mobil_' gFile{gg} '.png'])
end