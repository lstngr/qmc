nbins = [10 50 100 500 1000];
filelist = "bin" + strtrim(cellstr(num2str(nbins'))') + ".out";
data = cell(1,length(filelist));

for idx=1:length(filelist)   
    formatSpec = '%f%f';
    fileID = fopen(filelist(idx),'r');
    data{idx} = textscan(fileID,formatSpec,'CommentStyle','#','Delimiter',';');
end