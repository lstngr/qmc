pvals = ceil(linspace(10,100,10));
filelist = run_scan('InputConfig/default.in','beta',pvals);
data = cell(1,length(filelist));

for idx=1:length(filelist)   
    formatSpec = '%f%f';
    fileID = fopen(filelist{idx},'r');
    data{idx} = textscan(fileID,formatSpec,'CommentStyle','#','Delimiter',';');
end
true_energy=-5.65105; % 8 spins (beta20)