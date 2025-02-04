function [outMat,count] = groupbypermute(inMat)

%Empty output matrix
outMat = NaN(size(inMat,1),1);

%Number of dimensions to group by
nDims = size(inMat,2);

%Combo is the set of all possible combinations
command = 'combo = allcomb(';
command2 = 'isIn = find(';
for nn = 1:nDims
    %Find the unique values in each dimension; assign to a variable
    eval(['valList' num2str(nn) '= unique(inMat(~isnan(inMat(:,nn)),nn));'])
    
    %Build the analysis commands
    if nn == 1
        command = [command 'valList' num2str(nn)];
        command2 = [command2 'inMat(:,1) == combo(ii,1)'];
    else
        command = [command ', valList' num2str(nn)];
         command2 = [command2 '& inMat(:,' num2str(nn) ') == combo(ii,' ...
             num2str(nn) ')'];
    end
end
command = [command ');'];
command2 = [command2 ');'];

%This will create the full combination list
eval(command)

%Go through each of the combinations and determine which group everything
%is in
for ii = 1:size(combo,1)
    eval(command2)
    outMat(isIn) = ii;
end

%Counts
[unis] = unique(outMat(~isnan(outMat)));

count = NaN(length(unis),2);
for ii = 1:length(unis)
    count(ii,1) = unis(ii);
    count(ii,2) = length(find(outMat == unis(ii)));
end