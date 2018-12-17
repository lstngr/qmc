% Runs a parameter scan within the range pvals
pvals = 2.^linspace(4,9,6);
filelist = run_scan('InputConfig/default.in','beta',pvals);
data = cell(1,length(filelist));

% Imports single measurements in a container (data)
for idx=1:length(filelist)
    data{idx} = import_single(filelist{idx},3);
end