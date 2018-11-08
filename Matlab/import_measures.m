filelist = "num_sweep=" + strtrim(cellstr(num2str(pvals','%.4f'))') + ".out";
data = cell(1,length(filelist));

for idx=1:length(filelist)   
    formatSpec = '%f%f';
    fileID = fopen(filelist{idx},'r');
    data{idx} = textscan(fileID,formatSpec,'CommentStyle','#','Delimiter',';');
end
true_energy=-5.65105; % 8 spins (beta20)