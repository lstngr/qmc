% Runs a parameter scan within the range pvals
% and imports the resulting data in the workspace
pvals = 2.^linspace(4,9,6);
filelist = run_scan('InputConfig/default.in','beta',pvals);
data = cell(1,length(filelist));

% Imports single measurements in a container (data)
% This container can be passed to browse_import.
for idx=1:length(filelist)
    data{idx} = import_single(filelist{idx},3);
end