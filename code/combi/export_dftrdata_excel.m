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

% check for missing pids in filenames
for i=1:56
    pid_prefix = ['PP' num2str(i,'%02d')];
    midx = find(contains(filenames,pid_prefix));
    if isempty(midx)
%         fprintf('Empty pid: %d \n',i);
    end
end

% check treatment labels
nb_pids = length(placebo.labels.pid);
for i=1:nb_pids
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


%% reformat to table, save as excel file

varnames = [{'condition'},attachment.names,log.ftr_names];
varnames{2} = 'SubjID';
data = cell(56*2,length(varnames));

for i=1:56
    pid_prefix  = ['PP' num2str(i,'%02d')];
    pattern_RS  = pid_prefix + "_RS_";
    pattern_EYE = pid_prefix + "_EYE_";
    
    midx_RS  = find(contains(filenames,pattern_RS));
    midx_EYE = find(contains(filenames,pattern_EYE));

    % RS data: fill first [nb_pids] rows with [diff(post,pre)]
    data = get_firstcols(data,i,attachment,'RS');
    data = get_dftrs(data,i,midx_RS,log,'RS',nb_pids);
    
    % EYE data: fill next [nb_pids] rows with [diff(post,pre)]
    data = get_firstcols(data,i,attachment,'EYE');
    data = get_dftrs(data,i,midx_EYE,log,'EYE',nb_pids);
end

TableFtr = cell2table(data,'VariableNames',varnames);
filename = [root 'combi\dftr_data.xlsx'];

writetable(TableFtr,filename,'Sheet',1);
test = readtable(filename);


%% Local functions

function data = get_firstcols(data,i,attachment,condition)
% fill the first columns with general information
if strcmp(condition,'EYE'), j = i+56; else j=i; end
data{j,1} = condition;
for k=1:5
    data{j,k+1} = attachment.scores(i,k);
end
end

function data = get_dftrs(data,i,midx,log,condition,nb_pids)
% collect difference of feature values pre to post

if strcmp(condition,'EYE')
    i = i+nb_pids;
end

if length(midx)==2
    for j=1:length(log.ftr_names)
        val_pst = log.ftr_values(midx(1),j);
        val_pre = log.ftr_values(midx(2),j);
        data{i,j+6} = val_pst - val_pre;
    end
elseif isempty(midx)
    for j=1:length(log.ftr_names)
        data{i,j+6} = nan;
    end
else
    disp('Warning - unexpected indices found')
end
end
