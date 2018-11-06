function filelist = run_scan(conffile,param,values)
    
    param = string(param);
    filelist = cell(1,length(values));
    
    for idx=1:length(values)
        v = values(idx);
        if ~isstring(v)
            v = num2str(v,'%.4f');
        end
        filelist{idx} = param+"="+v+".out";
        % Modify configuration file
        sed_opr = "sed -ri 's/^("+param+"=)(.*)$/\1'"+v+"'/g' "+conffile;
        sed_out = "sed -ri 's/^(output_dest=)(.*)$/\1'"+filelist{idx}+"'/g' "+conffile;
        system(sed_opr);
        system(sed_out);
        
        % Run QMC
        sys_runcommand = "./qmc "+conffile;
        system(sys_runcommand);
    end
    
end