function [dftr,rownames,colnames] = dftr_mat(log,filenames,placebo)
% implementation based on export_ftrdata_excel.m
% generate treatment effect matrix from logdata

nb_pids = length(placebo.labels.pid);
dftr = nan(nb_pids*2,length(log.ftr_names));

colnames = log.ftr_names;
rownames = cell(1,nb_pids*2);

for i=1:nb_pids
    pid_prefix  = ['PP' num2str(i,'%02d')];
    pattern_RS  = pid_prefix + "_RS_";
    pattern_EYE = pid_prefix + "_EYE_";
    
    midx_RS  = find(contains(filenames,pattern_RS));
    midx_EYE = find(contains(filenames,pattern_EYE));

    % RS data: fill first nb_pids rows with [diff(post,pre)]
    [dftr,rownames] = get_dftrs(dftr,rownames,i,midx_RS,log,nb_pids,'RS');
    
    % EYE data: fill next nb_pids rows with [diff(post,pre)]
    [dftr,rownames] = get_dftrs(dftr,rownames,i,midx_EYE,log,nb_pids,'EYE');
end
end


%% Local function

function [dftr,rownames] = get_dftrs(dftr,rownames,i,midx,log,nb_pids,condition)
% collect difference of feature values pre to post

    rowname = ['PP' num2str(i,'%02d') '_' condition];
    if strcmp(condition,'EYE')
        i = i+nb_pids;
    end
    rownames{i} = rowname;

    if length(midx)==2
        for j=1:length(log.ftr_names)
            val_pst = log.ftr_values(midx(1),j);
            val_pre = log.ftr_values(midx(2),j);
            dftr(i,j) = val_pst - val_pre;
        end
    elseif isempty(midx)
        for j=1:length(log.ftr_names)
            dftr(i,j) = nan;
        end
    else
        disp('Warning - unexpected indices found')
    end
end