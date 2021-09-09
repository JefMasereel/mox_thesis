close all
clear

%% Feature analysis - apply tests and build result reports
% using TFDATA from build_tfdata.m (transformations applied)

% Stress response - finding reliable features for further analysis
% expanded to 4 tests of equal sample size for fairer comparison

% 1a. pre,oxt,rs/eye  |
% 1b. pre,plc,rs/eye  |
% 2. post,oxt,rs/eye  |> pval, pFDR -> list sign. ftrs ranked by pFDR
% 3. post,plc,rs/eye  |


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
orglog      = load([datadir 'logdata_all.mat']); % original ftr data
placebo     = load([datadir 'placebo.mat']);
attachment  = load([datadir 'attachment.mat']).attachment;
filenames   = load([datadir 'filenames.mat']).filenames;

% load prepped data from build_tfdata.m
SOURCE = load([dir_intermediates 'normtf_data.mat']);
tf_info = SOURCE.transformations;
ftr_ref = SOURCE.ftr_ref;
DATA = SOURCE.TFDATA;

% build data handling structures & treatment effect data (dftr)
% write new copy of log.mat to use prepped DATA iso. raw features
log = orglog; log.ftr_values = DATA;
hdl = build_handlers(log,placebo,attachment,filenames);
pids = hdl.pids;
idx  = hdl.idx;


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


%% General specs

% inspection: use as loaded from SOURCE
inspect.normality = SOURCE.normality;
inspect.vartestn = SOURCE.vartestn;

% main test specs
test.data   = DATA;
test.paired = true;
test.method = 'signrank';
test.alpha = .05;

% toggle output
verbose = true;  % print assumption check results
plotFDR = false; % individual ranked pFDR plots
FDR_all = true;  % grouped figure with ranked results


%% Prepare table construction

% working with ftr values loaded from DATA
% apply 3 tests, get pvals & pFDR from each
% sort rows by pFDR ranking, only include if p<alpha

nb_features = size(test.data,2);

% table layout = {p1 p2 p3 pFDR1 pFDR2 pFDR3}
tbl.pval = cell(4,1);
tbl.pFDR = cell(4,1);
tbl.varnames = {'p1','p2','p3','p4','q1','q2','q3','q4'};
tbl.rownames = log.ftr_names;



%% Stress response at baseline
% pre-treatment, oxt only
% does the feature set capture a significant stress response?

test.idx = [(idx.PRE & idx.oxytocin & idx.RS),...
            (idx.PRE & idx.oxytocin & idx.EYE)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);

% estimate pFDR across features [Storey2002]
% or use 'BHFDR',true for BH method
[test.FDR,test.q] = mafdr(test.p);

% review results
if plotFDR
    % large ranked barplot
    plot_rankedFDR(test,log);
end

% save relevant results to Table
tbl.pval{1} = test.p;
tbl.qval{1} = test.q;
tbl.pFDR{1} = test.FDR;


%% Stress response at baseline
% pre-treatment, plc only
% does the feature set capture a significant stress response?

test.idx = [(idx.PRE & idx.placebo & idx.RS),...
            (idx.PRE & idx.placebo & idx.EYE)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);

% estimate pFDR across features [Storey2002]
% or use 'BHFDR',true for BH method
[test.FDR,test.q] = mafdr(test.p);

% review results
if plotFDR
    % large ranked barplot
    plot_rankedFDR(test,log);
end

% save relevant results to Table
tbl.pval{2} = test.p;
tbl.qval{2} = test.q;
tbl.pFDR{2} = test.FDR;


%% Stress response after oxytocin treatment
% does the feature set capture a significant stress response?

test.idx = [(idx.POST & idx.oxytocin & idx.RS),...
            (idx.POST & idx.oxytocin & idx.EYE)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);

% estimate pFDR across features [Storey2002]
% or use 'BHFDR',true for BH method
[test.FDR,test.q] = mafdr(test.p);

% review results
if plotFDR
    % large ranked barplot
    plot_rankedFDR(test,log);
end

% save relevant results to Table
tbl.pval{3} = test.p;
tbl.qval{3} = test.q;
tbl.pFDR{3} = test.FDR;


%% Stress response after placebo treatment
% does the feature set capture a significant stress response?

test.idx = [(idx.POST & idx.placebo & idx.RS),...
            (idx.POST & idx.placebo & idx.EYE)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);

% estimate pFDR across features [Storey2002]
% or use 'BHFDR',true for BH method
[test.FDR,test.q] = mafdr(test.p);

% review results
if plotFDR
    % large ranked barplot
    plot_rankedFDR(test,log);
end

% save relevant results to Table
tbl.pval{4} = test.p;
tbl.qval{4} = test.q;
tbl.pFDR{4} = test.FDR;
% tbl.sampsize{3} = 


%% Review results in matlab table

tbl_mat = table(tbl.pval{1},tbl.pval{2},tbl.pval{3},tbl.pval{4},...
                tbl.qval{1},tbl.qval{2},tbl.qval{3},tbl.qval{4},...
                'VariableNames',tbl.varnames,...
                'RowNames',tbl.rownames);

% sort by first col: pval(pre)
[tbl_mat_sorted,I] = sortrows(tbl_mat); 

% only include significant rejections
filter_sorted = (tbl_mat_sorted.p1 < test.alpha) | ...
                (tbl_mat_sorted.p2 < test.alpha) | ...
                (tbl_mat_sorted.p3 < test.alpha) | ...
                (tbl_mat_sorted.p4 < test.alpha) | ...
                (tbl_mat_sorted.q1 < test.alpha) | ...
                (tbl_mat_sorted.q2 < test.alpha) | ...
                (tbl_mat_sorted.q3 < test.alpha) | ...
                (tbl_mat_sorted.q4 < test.alpha);
tbl_tex = tbl_mat_sorted(filter_sorted,:);

if verbose
    disp(tbl_tex)
end

% ftrIDs sorted & filtered for compact view
tex.fID_sort = ftr_ref.ID_str(I);
tex.fID_filt = ftr_ref.ID_str(filter_sorted);

% printable header names
tex.rows = strrep(tbl_tex.Properties.RowNames,'_',' ');
tex.cols = tbl_tex.Properties.VariableNames;
tex.data = tbl_tex.Variables;


%% build latex table

fprintf('Latex table\n\n');

input.tableCaption = 'Significant stress response features';
input.tableLabel = 'stress_response';
input.tableBorders = 0;

input.data = tex.data;
input.tableRowLabels = tex.rows; % ftr_ref.ID_str(filter);
input.tableColLabels = tex.cols;

input.tableColumnAlignment = 'r';
input.dataFormat = {'%.4f',8};

latex = latexTable(input);


%% Alternative summarized format
% compressed format compared to other scripts due to size of results

alt.rownames = strrep(pad(log.ftr_names),'_',' ');
alt.varnames.headers = {'feature','sample','pval','qval'};
alt.varnames.test = {'pre_oxt','pre_plc','post_oxt','post_plc'};

for i = 1:4         % sample / test number
    
%     fprintf('\n');
    fprintf('\\hline \n');

    [tmp.pval,I] = sort(tbl.pval{i});
    filter = find(tmp.pval < test.alpha);

    tmp.qval = tbl.qval{i}(I);
    tmp.ftr  = alt.rownames(I);
    tmp.test = alt.varnames.test{i};

    for j = 1:length(filter)

        fprintf('%s & %s & %.4f & %.4f %s \n',...
            tmp.ftr{j},tmp.test,...
            tmp.pval(j),tmp.qval(j),'\\');
    end
end

fprintf('\n\n');


%% Build summary figure

txt = {'pre, oxyotocin','pre, placebo','post, oxytocin','post, placebo'};
ftrnames = strrep(log.ftr_names,'_',' ');

if FDR_all
%     figure
%     title('Feature results in stress response tests')
%     for i=1:3
%         subplot(3,1,i)
%         subplot_rankedFDR(tbl.pval{i},tbl.pFDR{i},test.alpha,ftr_ref.ID_str,txt{i},'all');
%     end
    
    figure
    title('Feature results in stress response tests')
    for i=1:4
        subplot(4,1,i)
        subplot_ranked_pvals(tbl.pval{i},tbl.qval{i},test.alpha,...
                             ftrnames,txt{i},'sign_p'); % ftr_ref.ID_str
    end
    subplot(4,1,1)
    legend('p-values','q-values','alpha')

%     figure
%     title('Overview of significant features in stress response tests')
%     for i=1:3
%         subplot(3,1,i)
%         subplot_rankedFDR(tbl.pval{i},tbl.pFDR{i},test.alpha,ftr_ref.ID_str,txt{i},41);
%     end
end
