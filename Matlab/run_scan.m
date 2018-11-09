function filelist = run_scan(conffile,param,values)
    
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