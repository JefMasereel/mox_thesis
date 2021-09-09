function BCSR = build_bcsr(log,placebo)
% cfr. build_dftr.m, but for baseline-corrected stress response (bcsr)
% indices ridx.OT/PL/PRE/POST constructed in build_handlers.m

    nb_pids = length(placebo.labels.pid);
    bcsr = nan(nb_pids*2,length(log.ftr_names));
    
    filenames = log.filenames;    % row labels, input  (ftr)
    rownames = cell(1,nb_pids*2); % row labels, output (bcsr)
    colnames = log.ftr_names;     % col labels: features

    for i=1:nb_pids
        pid_prefix  = ['PP' num2str(i,'%02d')];
        idx.pid  = startsWith(filenames,pid_prefix);
        idx.pre  = contains(filenames,'_PRE');
        idx.post = contains(filenames,'_POST');

        midx.PRE  = find(idx.pid & idx.pre);
        midx.POST = find(idx.pid & idx.post);
        
        % PRE data: fill first nb_pids rows with [diff(eye,rst)]
        [bcsr,rownames] = get_bcsr(bcsr,rownames,i,midx.PRE,log,nb_pids,'PRE');

        % POST data: fill next nb_pids rows with [diff(eye,rst)]
        [bcsr,rownames] = get_bcsr(bcsr,rownames,i,midx.POST,log,nb_pids,'POST');
    end

    % return baseline-corrected stress response (bc_eye) as structure
    BCSR.data = bcsr;
    BCSR.rownames = rownames;
    BCSR.colnames = colnames;
end


%% Local function

function [bcsr,rownames] = get_bcsr(bcsr,rownames,i,midx,log,nb_pids,state)
% collect difference of feature values rs/eye

    rowname = ['PP' num2str(i,'%02d') '_' state];
    if strcmp(state,'POST')
        i = i+nb_pids;
    end
    rownames{i} = rowname;

    if length(midx)==2
        for j=1:length(log.ftr_names)
            % assuming filenames are ordered!
            % midx(1) = eye, midx(2) = rest
            
            % debug
            assert(contains(log.filenames(midx(1)),'_EYE_'));
            assert(contains(log.filenames(midx(2)),'_RS_'));
            
            val_eye = log.ftr_values(midx(1),j);
            val_rst = log.ftr_values(midx(2),j);
            bcsr(i,j) = val_eye - val_rst;
        end
    elseif isempty(midx)
        fprintf('Warning - no data available for %s\n',rowname);
        bcsr(i,:) = nan;
    else
        fprintf('Warning - incomplete pair for %s\n',log.filenames{midx});
        bcsr(i,:) = nan;
    end
end