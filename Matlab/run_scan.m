function filelist = run_scan(conffile,param,values)
% Runs a multi-process parameter scan, the name of the parameter
% is provided in 'param', its values in 'values'. The remaining QMC
% parameters are fetched from the configuration file 'conffile'.
% The function generates temporary configuration file in the root folder,
% and stores simulation outputs in './out'.
% In the event of existing output files, the SSE QMC data is appended and
% split by the init line of the program.

% NOTE: This script calls a parfor (parallel-for) to run multiple
% simulations at the same time. In the event that the 'Parallel Computing
% Toolbox' is unavalaible on the client machine, a usual for-loop will be
% substituted.

% NOTE: One can call run_scan with param=num_sweep, but the doing so just 
% generates multiple output files with the same data at the head. The
% process is equivalent to running one simulation for the longest num_sweep
% desired and truncating the output file at the desired locations. run_scan
% may still be used for convenience if the simulation time stays reasonably
% low.

param = string(param);
filelist = cell(1,length(values));
parfor_conf = filelist;
if ~exist('out','dir')
    !mkdir out/
end

parfor idx=1:length(values)
    v = values(idx);
    if ~isstring(v)
        v = num2str(v,'%.4f');
    end
    filelist{idx} =param+"="+v+".out";
    % Modify configuration file
    parfor_conf{idx} = "parfor"+num2str(idx)+".in";
    system("cp "+conffile+" "+parfor_conf{idx});
    sed_opr = "sed -ri 's/^("+param+"=)(.*)$/\1'"+v+"'/g' "+parfor_conf{idx};
    sed_out = "sed -ri 's/^(output_dest=)(.*)$/\1'"+"'out\/'"+filelist{idx}+"'/g' "+parfor_conf{idx};
    system(sed_opr);
    system(sed_out);
    
    % Run QMC
    sys_runcommand = "./qmc "+parfor_conf{idx};
    disp(sys_runcommand)
    system(sys_runcommand);
    system("rm "+parfor_conf{idx})
    
    % Update result path
    filelist{idx} = "out/"+filelist{idx};
end

end