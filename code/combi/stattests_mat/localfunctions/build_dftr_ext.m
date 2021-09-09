function DFTR = build_dftr_ext(log,attachment,filenames,placebo)
% implementation based on export_ftrdftr_excel.m and dftr_mat.m
% generate treatment effect matrix from logdftr
% WITH extended columns to include attachment scores.

    nb_pids = length(placebo.labels.pid);
    dftr = cell(nb_pids*2,length(log.ftr_names));

    colnames = [{'condition'},attachment.names,log.ftr_names];
    rownames = cell(1,nb_pids*2);

    for i=1:nb_pids
        pid_prefix  = ['PP' num2str(i,'%02d')];
        pattern_RS  = pid_prefix + "_RS_";
        pattern_EYE = pid_prefix + "_EYE_";

        midx_RS  = find(contains(filenames,pattern_RS));
        midx_EYE = find(contains(filenames,pattern_EYE));

        % RS dftr: fill first nb_pids rows with [diff(post,pre)]
        dftr = get_firstcols(dftr,i,attachment,nb_pids,'RS');
        [dftr,rownames] = get_dftrs(dftr,rownames,i,midx_RS,log,nb_pids,'RS');

        % EYE dftr: fill next nb_pids rows with [diff(post,pre)]
        dftr = get_firstcols(dftr,i,attachment,nb_pids,'EYE');
        [dftr,rownames] = get_dftrs(dftr,rownames,i,midx_EYE,log,nb_pids,'EYE');
    end


    % return dftr as structure
    DFTR.data = dftr;
    DFTR.rownames = rownames;
    DFTR.colnames = colnames;
end


%% Local functions

function dftr = get_firstcols(dftr,i,attachment,nb_pids,condition)
% fill the first columns with general information
if strcmp(condition,'EYE'), j = i+nb_pids; else j=i; end
dftr{j,1} = condition;
for k=1:5
    dftr{j,k+1} = attachment.scores(i,k);
end
end

function [dftr,rownames] = get_dftrs(dftr,rownames,i,midx,log,nb_pids,condition)
% collect difference of feature values pre to post

    rowname = ['PP' num2str(i,'%02d') '_' condition];
    if strcmp(condition,'EYE')
        i = i+nb_pids;
    end
    rownames{i} = rowname;

    if length(midx)==2
        for j=1:length(log.ftr_names)
            % assuming filenames are ordered!
            % midx(1) = post, midx(2) = pre
            val_pst = log.ftr_values(midx(1),j);
            val_pre = log.ftr_values(midx(2),j);
            dftr{i,j+6} = val_pst - val_pre;
        end
    elseif isempty(midx)
        for j=1:length(log.ftr_names)
            dftr{i,j+6} = nan;
        end
    else
        disp('Warning - unexpected indices found')
    end
end