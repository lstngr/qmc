pvals = ceil(logspace(3,7,9));
filelist = run_scan('InputConfig/default.in','num_sweep',pvals);
data = cell(1,length(filelist));

for idx=1:length(filelist)   
    formatSpec = '%f%f';
    fileID = fopen(filelist{idx},'r');
    data{idx} = textscan(fileID,formatSpec,'CommentStyle','#','Delimiter',';');
end
true_energy=-5.65109; % 8 spins (beta60)