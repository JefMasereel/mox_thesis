close all
clear

%% Feature analysis - preparation of extracted features



%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\';
datadir = [root 'data\'];

% save/load intermediate results in data folder
dir_intermediates = [datadir 'combi\intermediates_analysis\'];

% utility functions
% localfunctions self-defined within repo
% https://github.com/bastibe/Violinplot-Matlab
% https://github.com/eliduenisch/latexTable
addpath([root 'code\combi\stattests_mat\localfunctions\']);
addpath('C:\Users\jefma\repos\Violinplot-Matlab');
addpath('C:\Users\jefma\repos\latexTable');


%% Load data from previous steps and prepare for testing formats

% load relevant source data
log         = load([datadir 'logdata_all.mat']); % original ftr data
placebo     = load([datadir 'placebo.mat']);
attachment  = load([datadir 'attachment.mat']).attachment;
filenames   = load([datadir 'adjusted_filenames.mat']).filenames;

% build data handling structures & treatment effect data (dftr)
hdl = build_handlers(log,placebo,attachment,filenames);
pids = hdl.pids;
idx = hdl.idx;

% define reference ID numbers for feature names
ftr_ref.ID_nb  = 1:length(log.ftr_names);
ftr_ref.ID_str = compose('%2d',ftr_ref.ID_nb);
ftr_ref.names  = log.ftr_names;


%% Test implementations
% general assumption checks & normality transformations
% apply transformations to a ftr if normality is reached
% apply zscore standardization to all features
% save resulting data to normtf_data.mat


%% Notes on test and sample specification in this script
% applied some OOP-ish variable construction to clean up testing process

% demo test specification
% test.paired = true;
% test.method = 'ttest';
% test.data   = log.ftr_values;
% test.idx    = [(idx.PRE & idx.RS),... % sample A
%                (idx.PRE & idx.EYE)];  % sample B

% note: test.idx should contain logical values!
% using ones and zeros breaks the indexing..

% test.method =
% ttest           parametric      paired      difference
% signrank        non-parametric  paired      difference
% ttest2          parametric      non-paired  difference
% ranksum         non-parametric  non-paired  difference
% pearson_corr    parametric      ?           correlation
% spearman_corr   non-parametric  ?           correlation


%% General assumption checks and normality transformation attempts
% only apply changes to full dataset! avoid sample differentiation

% full dataset as single sample
test.paired = false;
test.method = NaN;
test.data = log.ftr_values;
test.idx = true(size(idx.RS));

% specify assumption checks
inspect.vartestn.method     = 'Bartlett'; % less strict: OBrien, Levene..
inspect.vartestn.alpha      = .05;
inspect.normality.method    = 'jbtest';   % unknown mean & variance!
inspect.normality.alpha     = .05;

% save test specs to check the same assumptions on smaller samples
% do NOT apply any transformations on subsamples! avoid artf separation
ORG_INSPECT   = inspect;


% specify normality transformation methods
% notes: 
% - transformations are applied in fixed order!!
% - add intra-subject standardization option? 
%   would be better to move this to ftr extraction step if included
tf.boxcox.apply          = true;
tf.rmv_outliers.apply    = true;
tf.rmv_outliers.method   = 'mean';
tf.zscore.apply          = false;       % currently applied after!


% suppress auto-report prints in this section, not useful. 
verbose=false; 

% apply
test          = build_samples(test,pids,filenames);
inspect       = inspect_distribution(test,inspect,verbose);
[tf_test,tf]  = apply_transformations(test,tf,filenames);
tf_inspect    = inspect_distribution(tf_test,inspect,verbose);

% vartestn assumptions not relevant for single sample! ignore auto-report
ratio_before = 100*sum(inspect.normality.h)/test.nb_features;
ratio_after  = 100*sum(tf_inspect.normality.h)/tf_test.nb_features;
fprintf("Normality tests:\n%.2f%% rejected before tf\n%.2f%% rejected after tf\n",...
            ratio_before,ratio_after);


%% Review test results in table format

% data missingness
% test.nb_nans : missing values directly from feature extraction
nb_rmv_outliers = tf_inspect.nb_nans-inspect.nb_nans;

RowNames = log.ftr_names;
VarNames = {'original','transformed','invalid','outliers','lambda','offset'};
tf_overview = table(inspect.normality.p, tf_inspect.normality.p, ...
                    inspect.nb_nans,nb_rmv_outliers,...
                    tf.boxcox.lambda, tf.boxcox.offset, ...
                    'RowNames',RowNames,'VariableNames',VarNames);

% disp(tf_overview) % full table
tmp_filter = tf_overview.transformed>inspect.normality.alpha;
disp(tf_overview(tmp_filter,:))
fprintf('%d features fit normal distribution after transformations\n',sum(tmp_filter));

% printable header names
tex.rows = strrep(tf_overview.Properties.RowNames,'_',' ');


%% write table to latex

fprintf('Latex table\n\n');

input.tableCaption = 'Overview of applied transformations';
input.tableLabel = 'tf_overview';
input.tableBorders = 0;

input.tableColumnAlignment = 'c';
input.data = tf_overview(tmp_filter,:).Variables; % array
input.dataFormat = {'%.4f',2,'%d',2,'%.4f',1,'%d',1};
input.tableRowLabels = tex.rows(tmp_filter); %ftr_ref.ID_str(tmp_filter);
input.tableColLabels = tf_overview.Properties.VariableNames;

latex = latexTable(input);


%% push useful transformations back to matrix format for further tests

% assumes test.sample contains only one sample!
assert(tf_test.nb_samples==1);
tf_data = nan(size(tf_test.data));
for ftr = 1:tf_test.nb_features
    if tf_inspect.normality.h(ftr)==0
        tf_data(:,ftr) = tf_test.sample{ftr}(:,1); % use transformed data
    elseif tf_inspect.normality.h(ftr)==1
        tf_data(:,ftr) = test.sample{ftr}(:,1);    % use original data
    end
end

% apply zscore to all features, ignoring NaN values
tf_data = custom_zscore(tf_data);
tf.zscore.apply = true;


%% save transformed data for later use and logging purposes

% relevant variables to save
normality = rmfield(tf_inspect.normality,{'p'});
vartestn = rmfield(tf_inspect.vartestn,{'p'});
transformations = tf;
TFDATA = tf_data;

clear input
prompt = '\nSave these transformed data for further analysis? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'N';
end
if strcmp(str,'Y')
    savepath = [dir_intermediates 'normtf_data.mat'];
    save(savepath, 'normality','vartestn','transformations','TFDATA','ftr_ref');
    fprintf('Results saved to %s\n',savepath);
else
    fprintf('Results not saved.\n');
end


%% review missing data

pattern = 'RSP_RR_MN';    % search pattern for feature names or groups
threshold = 0;              % minimum nb of nans
nb_elm = numel(log.ftr_values(:,find(contains(log.ftr_names,pattern))));
ratio = sum(isnan(log.ftr_values(:,find(contains(log.ftr_names,pattern)))),'all') / nb_elm*100;

% invalid from logfile
fprintf('%s: %.2f%% of datapoints is empty (nan)\n',pattern,ratio);
disp(log.filenames(find(sum(isnan(log.ftr_values(:,find(contains(log.ftr_names,pattern)))),2)>threshold))) 

% mids removed as outliers
disp(tf.rmv_outliers.mids{contains(log.ftr_names,pattern)})

% nb_elm = numel(log.ftr_values(:,find(contains(log.ftr_names,'CDA'))));
% ratio = sum(isnan(log.ftr_values(:,find(contains(log.ftr_names,'CDA')))),'all') / nb_elm*100;
% fprintf('CDA: %.2f%% of datapoints is empty (nan)\n',ratio);
% disp(log.filenames(find(sum(isnan(log.ftr_values(:,find(contains(log.ftr_names,'DDA')))),2)>0)))
% 
% nb_elm = numel(log.ftr_values(:,find(contains(log.ftr_names,'DDA'))));
% ratio = sum(isnan(log.ftr_values(:,find(contains(log.ftr_names,'DDA')))),'all') / nb_elm*100;
% fprintf('DDA: %.2f%% of datapoints is empty (nan)\n',ratio);
% disp(log.filenames(find(sum(isnan(log.ftr_values(:,find(contains(log.ftr_names,'DDA')))),2)>0)))
% 
% nb_elm = numel(log.ftr_values(:,find(contains(log.ftr_names,'RSP'))));
% ratio = sum(isnan(log.ftr_values(:,find(contains(log.ftr_names,'RSP')))),'all') / nb_elm*100;
% fprintf('RSP: %.2f%% of datapoints is empty (nan)\n',ratio);
% disp(log.filenames(find(sum(isnan(log.ftr_values(:,find(contains(log.ftr_names,'RSP')))),2)>0)))


%% old versions
% % % print select features (manual selection)
% % fprintf('%d empty mids in PRV features: \n',sum(isnan(log.ftr_values(:,1))));
% % disp(log.filenames(find(isnan(log.ftr_values(:,1)))));
% % 
% % fprintf('%d empty mids in CDA IRI features: \n',sum(isnan(log.ftr_values(:,17))));
% % disp(log.filenames(find(isnan(log.ftr_values(:,17)))));
% % 
% % fprintf('%d empty mids in DDA IRI features: \n',sum(isnan(log.ftr_values(:,28))));
% % disp(log.filenames(find(isnan(log.ftr_values(:,28)))));
% 
% % % iterate over features
% % for ftr = 1:size(log.ftr_values,2)
% %     disp(ftr);disp(log.filenames(find(isnan(log.ftr_values(:,ftr)))));
% % end
% 
% % % iterate over mids
% % for i = 1:size(log.ftr_values,1)
% %     nb_nans = sum(isnan(log.ftr_values(i,:)));
% %     if nb_nans > 5 % choose a threshold
% %         fprintf('mid %s contains %d empty features\n',log.filenames{i},nb_nans)
% %     end
% % end


%% review outliers if relevant

% % most outliers detected in DDA related features
% width_ftrs = find(contains(log.ftr_names,'width'));
% slope_ftrs = find(contains(log.ftr_names,'slope'));
% tf.rmv_outliers.mids(width_ftrs);
% tf.rmv_outliers.mids(slope_ftrs);
% temp = tf.rmv_outliers.mids(width_ftrs);
% table(temp{2},temp{3})


%% Local functions 

function [test,tf] = apply_transformations(test,tf,filenames)
% only to be used for complete dataset, not on separated subsamples!
    
    fprintf('applying sample transformations...\n');
    nb_samples  = test.nb_samples;
    nb_features = test.nb_features; 
    tf.boxcox.lambda = zeros(nb_features,nb_samples);
    tf.boxcox.offset = zeros(nb_features,nb_samples);
    tf.rmv_outliers.nb_removed = zeros(nb_features,1);
    tf.rmv_outliers.mids = cell(nb_features,1);

    for ftr = 1:nb_features
        
        % Box-Cox transformation
        if tf.boxcox.apply
            for i = 1:nb_samples
                tmp_org = test.sample{ftr}(:,i);
                tmp_offset = (ceil(-min(tmp_org)) + 1)*(min(tmp_org) <= 0);
                tmp_pos = tmp_org + tmp_offset;
                [tmp_bct,lambda] = boxcox(tmp_pos);
                test.sample{ftr}(:,i) = tmp_bct - tmp_offset;
                tf.boxcox.lambda(ftr,i) = lambda;
                tf.boxcox.offset(ftr,i) = tmp_offset;
            end
        end
        
        % outlier removal
        if tf.rmv_outliers.apply
            I = isoutlier(test.sample{ftr},tf.rmv_outliers.method);
            tf.rmv_outliers.nb_removed(ftr) = sum(I);
            tf.rmv_outliers.mids{ftr} = filenames(I);
            test.sample{ftr}(I) = nan;
        end
        
        % zscore standardization (intersubject)
        % removes mean from the samples!! loss of information!
        if tf.zscore.apply
            test.sample{ftr} = custom_zscore(test.sample{ftr});
        end
    end
end


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