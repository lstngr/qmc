nmeas  = zeros(1,length(data));
energy = nmeas; std_en = nmeas;

for i=1:length(data)
    
    pos_sweeps = data{i}{1};
    bins   = data{i}{2};
    
    outliers = isoutlier(bins);
    
%     plot(pos_sweeps(~outliers),bins(~outliers),'b-',...
%          pos_sweeps( outliers),bins( outliers),'rx')
%     xlabel 'nsweep'
%     ylabel 'Energy'
%     pause
%     
%     histfit(bins(~outliers),10,'norm')
%     pause
    [energy(i),std_en(i)] = normfit(bins(~outliers));
    nmeas(i) = length(bins(~outliers));
    
end