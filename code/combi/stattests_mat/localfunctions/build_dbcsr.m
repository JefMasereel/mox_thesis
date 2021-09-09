function DBCSR = build_dbcsr(bcsr,log,placebo)
% cfr. build_dftr.m, but for DBCSR = BCSR(post) - BCSR(pre)
% indices dridx.OT/PL constructed in build_handlers.m

    nb_pids = length(placebo.labels.pid);
    DBCSR.data = nan(nb_pids,length(log.ftr_names));
    
    filenames = bcsr.rownames;          % row labels, input  (bcsr)
    DBCSR.rownames = cell(1,nb_pids);   % row labels, output (dbcsr)
    DBCSR.colnames = log.ftr_names;     % col labels: features

    for i=1:nb_pids
        pid_prefix  = ['PP' num2str(i,'%02d')];
        DBCSR.rownames{i} = pid_prefix;
        
        idx.pid  = startsWith(filenames,pid_prefix);
        idx.pre  = contains(filenames,'_PRE');
        idx.post = contains(filenames,'_POST');

        midx.PRE  = find(idx.pid & idx.pre);
        midx.POST = find(idx.pid & idx.post);
        
        % DBCSR = BCSR(post) - BCSR(pre)
        if ~isempty(midx.PRE) && ~isempty(midx.PRE)
            for ftr = 1:length(log.ftr_names)
                val_pre = bcsr.data(midx.PRE,ftr);
                val_pst = bcsr.data(midx.POST,ftr);
                DBCSR.data(i,ftr) = val_pst - val_pre;
            end
        elseif isempty(midx.PRE) && isempty(midx.POST)
            fprintf('Warning - no data available for %s\n',pid_prefix);
            DBCSR.data(i,:) = nan;
        else
            fprintf('Warning - incomplete pair for %s\n',pid_prefix);
            DBCSR.data(i,:) = nan;
        end
    end
end