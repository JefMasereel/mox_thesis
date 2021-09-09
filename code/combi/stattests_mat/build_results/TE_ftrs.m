close all
clear

%% Feature analysis - apply tests and build result reports
% using TFDATA from build_tfdata.m (transformations applied)

% treatment effect - finding relevant features (pval,pFDR)
% 1. pre/post,rs,oxt
% 2. pre/post,rs,plc
% 3. pre/post,eye,oxt
% 4. pre/post,eye,plc
% _. pre/post,eye-rs,oxt
% _. pre/post,eye-rs,plc


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
idx = hdl.idx;


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
test.method = 'signrank'; %'ttest';
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
tbl.varnames = {'p1','q1',...
                'p2','q2',...
                'p3','q3',...
                'p4','q4'};
tbl.rownames = log.ftr_names;


%% Treatment effect - RS,OT
% does the feature set capture a significant treatment effect?

test.idx = [(idx.RS & idx.oxytocin & idx.PRE),...
            (idx.RS & idx.oxytocin & idx.POST)];

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


%% Treatment effect - RS,PL
% does the feature set capture a significant treatment effect?

test.idx = [(idx.RS & idx.placebo & idx.PRE),...
            (idx.RS & idx.placebo & idx.POST)];

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


%% Treatment effect - EYE,OT
% does the feature set capture a significant treatment effect?

test.idx = [(idx.EYE & idx.oxytocin & idx.PRE),...
            (idx.EYE & idx.oxytocin & idx.POST)];

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


%% Treatment effect - EYE,PL
% does the feature set capture a significant treatment effect?

test.idx = [(idx.EYE & idx.placebo & idx.PRE),...
            (idx.EYE & idx.placebo & idx.POST)];

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


%% Review results in matlab table

tbl_mat = table(tbl.pval{1},tbl.qval{1},...
                tbl.pval{2},tbl.qval{2},...
                tbl.pval{3},tbl.qval{3},...
                tbl.pval{4},tbl.qval{4},...
                'VariableNames',tbl.varnames,...
                'RowNames',tbl.rownames);

% tbl_mat = table(tbl.pval{1},tbl.pFDR{1},...
%                 tbl.pval{2},tbl.pFDR{2},...
%                 tbl.pval{3},tbl.pFDR{3},...
%                 tbl.pval{4},tbl.pFDR{4},...
%                 'VariableNames',tbl.varnames,...
%                 'RowNames',tbl.rownames);

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

input.tableCaption = 'Significant treatment effect features';
input.tableLabel = 'te_ftrs';
input.tableBorders = 0;

input.data = tex.data;
input.tableRowLabels = tex.rows; % tex.fID_filt; % 
input.tableColLabels = tex.cols;

input.tableColumnAlignment = 'c';
input.dataFormat = {'%.4f',8};

latex = latexTable(input);


%% Build summary figure

ftrnames = strrep(log.ftr_names,'_',' ');
txt = {'rest, oxt','rest, plc','stress, oxt','stress, plc'};
nb_tests = 4;

if FDR_all
%     figure
%     title('Feature results in stress response tests')
%     for i=1:nb_tests
%         subplot(nb_tests,1,i)
%         subplot_rankedFDR(tbl.pval{i},tbl.pFDR{i},test.alpha,ftr_ref.ID_str,txt{i},'all');
%     end
    
    figure
    title('Feature results in stress response tests')
    for i=1:nb_tests
        subplot(nb_tests,1,i)
        subplot_ranked_pvals(tbl.pval{i},tbl.qval{i},test.alpha,ftrnames,txt{i},'sign_p');
    end
    subplot(nb_tests,1,1)
    legend('p-values','q-values','alpha')

%     figure
%     title('Overview of significant features in stress response tests')
%     for i=1:nb_tests
%         subplot(nb_tests,1,i)
%         subplot_rankedFDR(tbl.pval{i},tbl.pFDR{i},test.alpha,ftr_ref.ID_str,txt{i},41);
%     end
end
