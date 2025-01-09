ncclear
mp.urlH = ...
    'dods://geoport.whoi.edu/thredds/dodsC/usgs/data2/sdalyander/SWAN_esp_tau/esp_agg.ncml';
nc = cfdataset(mp.urlH);

time = nc.time('time');
var = nc.variable('tauwc');

for in = 1:8
    start = (in-1)*1000+1;
    mEnd = in*1000;
    eval(['tauwc' num2str(in) '= double(var.data(' num2str(start) ':' ...
        num2str(mEnd) ',:,:));'])
end

tauwc9 = double(var.data(8001:8760,:,:));

tauwc = cat(1,tauwc1,tauwc2,tauwc3,tauwc4,tauwc5,tauwc6,tauwc7,tauwc8,tauwc9);
clear tauwc1 tauwc2 tauwc3 tauwc4 tauwc5 tauwc6 tauwc7 tauwc8 tauwc9
clear nc mEnd in var start

