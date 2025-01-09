%% Load data
lonNew = double(lon);
latNew = double(lat);
testNew = double(squeeze(m95(1,1,:,:)));
lonNew = lonNew(:);
latNew = latNew(:);
testNew = testNew(:);
lonNew(isnan(testNew)) = [];
latNew(isnan(testNew)) = [];
testNew(isnan(testNew)) = [];

DT = delaunay(lonNew,latNew);

%% Test TIN direct
fID = fopen('test.adf','w');
fwrite(fID,'TIN/n')
fwrite(fID,'BEGT/n')
fwrite(fID,'TNAM test/n')
fwrite(fID,'255 255 255/n')
fwrite(fID,['VERT ' num2str(length(testNew)) '\n'])

for ii = 1:length(testNew)
    fwrite(fID,[num2str(lonNew(ii)) ' ' num2str(latNew(ii)) ' ' ...
        num2str(testNew(ii)) ' 0\n'])
end
fwrite(fID,['TRI ' num2str(size(DT,1)) '\n'])
for ii = 1:size(DT,1)
    fwrite(fID,[num2str(DT(ii,1)) ' ' num2str(DT(ii,2)) ' ' num2str(DT(ii,3)) '\n'])
end
fwrite(fID,'ENDT')

%% Create xml
tic
x = lonNew;
y = latNew;
z = testNew;
xml_header = ['<?xml version="1.0"?>' ...
    '<LandXML version="1.1" date="2008-06-09" time="08:33:4"' ...
    'xmlns="http://www.landxml.org/schema/LandXML-1.1"' ...
    'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"' ...
    'xsi:schemaLocation="http://www.landxml.org/schema/LandXML-1.1 http://www.landxml.org/schema/LandXML-1.1/LandXML-1.1.xsd">' ...
    '<Units>' ...
    '<Metric areaUnit="squareMeter" linearUnit="meter" volumeUnit="cubicMeter"' ...
    'temperatureUnit="celsius" pressureUnit="milliBars"/>' ...
    '</Units>' ...
    '<Application name="FVCOM" version="10.2.0 (Build 400)" manufacturer="SMAST"' ...
    'manufacturerURL="http://fvcom.smast.umassd.edu/FVCOM/index.html"><Author createdBy="rsignell"/></Application>' ...    
    '<Surfaces>' ...
    '<Surface name="FVCOM_GOM2_Grid">'];   
f=fopen('test.xml', 'w');
fprintf(f,xml_header);
fprintf(f,['  <Definition surfType=\"TIN\" elevMin=\"' num2str(min(z(:))) '\" elevMax=\"' num2str(max(z(:))) '\">\n']);
fprintf(f, '    <Pnts>\n');
for i = 1:length(z)
    fprintf(f,['<P id=\"' num2str(i) '\"> ' num2str(x(i)) ' ' num2str(y(i)) ' ' num2str(z(i)) '</P>']);
end
fprintf(f,'   </Pnts>\n');
fprintf(f,'   <Faces>\n');
for i = 1:size(DT,1)
    fprintf(f,['<F> ' num2str(DT(i,1)) ' ' num2str(DT(i,2)) ' ' num2str(DT(i,3)) '</F>\n']);
end
fprintf(f,'   </Faces></Definition></Surface></Surfaces></LandXML>');
fclose(f);
toc
