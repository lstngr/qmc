% Reads a container (data), storing simulation data.
% data should be a cell containing outputs of import_single
nmeas  = zeros(1,length(data));
energy = nmeas; std_en = nmeas;
plotfigs = true; % Plot figures

for i=1:length(data)
    
    pos_sweeps = data{i}{1};
    bins   = data{i}{2}; % Energy bins
    
    % Default MATLAB routine, used to remove outliers
    % occuring in the first simulation sweeps. Occasionnaly
    % flags bins at larger times (see nmeas).
    outliers = isoutlier(bins);
    
    % Plot the energy vs the sweep, as well as the
    % distribution of the bins (fitted to a gaussian)
    if plotfigs
        plot(pos_sweeps(~outliers),bins(~outliers),'b-',...
            pos_sweeps( outliers),bins( outliers),'rx')
        xlabel 'nsweep'
        ylabel 'Energy'
        pause
        
        histfit(bins(~outliers),10,'norm')
        pause
    end

    [energy(i),std_en(i)] = normfit(bins(~outliers));
    % Nmeas is used to obtain the "correct" standard deviation of the sample:
    % std_en/sqrt(nmeas)
    nmeas(i) = length(bins(~outliers));
    
end