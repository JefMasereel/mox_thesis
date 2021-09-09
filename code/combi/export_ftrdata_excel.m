close all
clear

%% Feature analysis - export ftr data for Kaat & Nicky
% reformat logdata_all.mat and attachment.mat to ftr_data.xlsx
% as part of thesis by Jef Masereel


%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\data\';


tmp = load([root 'attachment.mat']);
attachment = tmp.attachment;

tmp = load([root 'adjusted_filenames.mat']);
filenames = tmp.filenames;
clear tmp

log = load([root 'logdata_all.mat']);
placebo = load([root 'placebo.mat']);


%% inspect data, clean up some labels

% check treatment labels
for i=1:56
    tmpa = attachment.scores(i,2);
    tmpb = placebo.labels.treatment{i};
%     assert((tmpa==1&strcmp(tmpb,'oxytocin'))|(tmpa==2&strcmp(tmpb,'placebo')));
%     % commented for demonstration with devalued data (empty pids)
end
disp('Treatment labels match')

% add 'PRV_' prefix to ftr values from classical PRV analysis
prv_ftrs = {'TOT','LF','HF','LFn','HFn','LH_ratio'};
for i=1:length(log.ftr_names)
    if any(strcmp(log.ftr_names{i},prv_ftrs))
        log.ftr_names{i} = ['PRV_' log.ftr_names{i}];
        disp(log.ftr_names{i})
    end
end

% check for missed ftrs without modality prefix
mod_prefix = {'RSPA','RSPC','PRV','CDA','DDA'};
for i=1:length(log.ftr_names)
    if any(strcmp(log.ftr_names{i},mod_prefix))
        disp(log.ftr_names{i})
    end
end
disp('all clear')


%% apply basic outlier removal (on individual feature columns)
% copied function from stattests_individual.m

% for each individual feature (univariate checks)
for i =1:length(log.ftr_names)
    [B,I] = rmoutliers(log.ftr_values(:,i),'quartiles');
    mask.outliers(:,i) = ~I;
    mask.rmv_counts(i) = sum(mask.outliers(:,i));
    mask.nan_counts(i) = sum(isnan(log.ftr_values(:,i)));
    
%     disp(log.ftr_names{i})
%     disp(sum(mask.outliers(:,i)))
end

% apply mask to overwrite outliers with NaN
newdata = nan(size(log.ftr_values));
for i =1:length(log.ftr_names)
    use_idx = find(mask.outliers(:,i)==1);
    newdata(use_idx,i) = log.ftr_values(use_idx,i);
end
log.ftr_values = newdata;

% quick review of removal counts:
disp(212 - sum(mask.outliers,1))


%% apply zscore standardization (comment if unwanted)
% copied function from stattests_individual.m

% for i = 1:length(placebo.labels.pid)
%     
%     pid_str = placebo.labels.pid{i};
%     pid_idx = find(contains(filenames,pid_str));
%     
%     % apply zscore standardization but ignore NaN values
%     log.ftr_values(pid_idx,:) = custom_zscore(log.ftr_values(pid_idx,:));
%     
% %     % use default zscore standardization (sensitive to NaNs)
% %     log.ftr_values(pid_idx,:) = zscore(log.ftr_values(pid_idx,:));
% end


%% reformat to table, save as excel file

varnames = [{'filename','session','condition'},attachment.names,log.ftr_names];
varnames{4} = 'SubjID';
varnames{5} = 'spray';
data = cell(length(log.filenames),length(varnames));

for i = 1:size(attachment.scores,1)
    pid_prefix = ['PP' num2str(i,'%02d')];
    midx = find(contains(filenames,pid_prefix));
    
    for j = midx'
        data{j,1} = filenames{j};
        data{j,2} = get_time(filenames{j});
        data{j,3} = get_cond(filenames{j});
        
        % export placebo labels as readable string (easier for plots)
        if attachment.scores(i,2)==1
            data{j,5} = 'oxytocin';
        elseif attachment.scores(i,2)==2
            data{j,5} = 'placebo';
        end
        
        for k=[1 3:5] % skip 2nd (placebo)
            data{j,k+3} = attachment.scores(i,k);
        end
        
        for k=1:length(log.ftr_names)
            data{j,k+8} = log.ftr_values(j,k);
        end
    end
end

TableFtr = cell2table(data,'VariableNames',varnames);

% save to file
filename = [root 'combi\ftr_data_raw_orm.xlsx']; % adjust filepath!!
writetable(TableFtr,filename,'Sheet',1);
test = readtable(filename);

disp('done')


%% Local functions

function Z = custom_zscore(A)
% apply zscore standardization to each column of matrix A
% BUT ignore any NaN values if present (iso. setting all to nan)

Z = nan(size(A));
for c = 1:size(A,2)
    column = A(:,c);
    column = column(~isnan(column));
    MN = mean(column);
    SD = std(column);
    for r = 1:size(A,1)
        if ~isnan(A(r,c))
            Z(r,c) = (A(r,c)-MN)/SD;
        end
    end
end
end

function strTime = get_time(filename)
if contains(filename,'_PRE')
    strTime = 'PRE';
elseif contains(filename,'_POST')
    strTime = 'POST';
end
end

function strCond = get_cond(filename)
if contains(filename,'_RS_')
    strCond = 'RS';
elseif contains(filename,'_EYE_')
    strCond = 'EYE';
end
end

