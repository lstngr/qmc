nbins = [1e3,5e3,1e4,5e4,1e5,5e5,1e6];
filelist = "num_sweep=" + strtrim(cellstr(num2str(nbins'))') + ".out";
data = cell(1,length(filelist));

for idx=1:length(filelist)   
    formatSpec = '%f%f';
    fileID = fopen(filelist(idx),'r');
    data{idx} = textscan(fileID,formatSpec,'CommentStyle','#','Delimiter',';');
end