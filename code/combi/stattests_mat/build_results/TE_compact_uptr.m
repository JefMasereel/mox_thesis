close all
clear

%% Feature analysis - apply tests and build result reports
% using TFDATA from build_tfdata.m (transformations applied)

% treatment effect - validating oxytocin wrt. placebo (pvals,pFDR)
% 1. dftr,rs,oxt/plc
% 2. dftr,eye,oxt/plc
% _. dftr,eye-rs,oxt/plc


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
didx = hdl.didx;
dftr = hdl.dftr;
ridx = hdl.ridx;
bcsr = hdl.bcsr;
dridx = hdl.dridx;
dbcsr = hdl.dbcsr;


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
test.paired = false;
test.method = 'ranksum'; %'ttest2';
test.alpha = .05;

% toggle output
verbose = true;  % print assumption check results
plotFDR = false; % individual ranked pFDR plots
FDR_all = true;  % grouped figure with ranked results


%% Prepare table construction

tbl.pval = cell(3,1);
tbl.pFDR = cell(3,1);
tbl.uptr = cell(3,1);
tbl.varnames = {'p1','q1','e1','p2','q2','e2','p3','q3','e3'};
tbl.rownames = log.ftr_names;


%% Treatment effect - dftr,RS,OT/PL
% do the features find significant differences in treatment effects?

test.data   = dftr.data;
test.idx = [(didx.RS & didx.oxytocin),...
            (didx.RS & didx.placebo)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,dftr.rownames);
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
tbl.uptr{1} = upwards_trend(test);


%% Treatment effect - dftr,EYE,OT/PL
% do the features find significant differences in treatment effects?

test.data   = dftr.data;
test.idx = [(didx.EYE & didx.oxytocin),...
            (didx.EYE & didx.placebo)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,dftr.rownames);
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
tbl.uptr{2} = upwards_trend(test);


%% Compare BCSR for treatment groups - DBCSR, oxt/plc

test.data   = dbcsr.data;
test.idx = [dridx.oxytocin,...
            dridx.placebo];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,dbcsr.rownames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);

% estimate pFDR across features [Storey2002]
% or use 'BHFDR',true for BH method
[test.FDR,test.q] = mafdr(test.p);

% review results
if plotFDR
    % large ranked barplot
    plot_rankedFDR(test,log);
    plot_rankedPvals(test,log);
end

% save relevant results to Table
tbl.pval{3} = test.p;
tbl.qval{3} = test.q;
tbl.pFDR{3} = test.FDR;
tbl.uptr{3} = upwards_trend(test);


%% Review results in matlab table

tbl_mat = table(tbl.pval{1},tbl.qval{1},tbl.uptr{1},...
                tbl.pval{2},tbl.qval{2},tbl.uptr{2},...
                tbl.pval{3},tbl.qval{3},tbl.uptr{3},...
                'VariableNames',tbl.varnames,...
                'RowNames',tbl.rownames);

% cleaner layout
[tbl_mat_sorted,I] = sortrows(tbl_mat); % sort by first col: pval(pre)
filter_sorted = (tbl_mat_sorted.p1 < test.alpha) | ...
                (tbl_mat_sorted.p2 < test.alpha) | ...
                (tbl_mat_sorted.p3 < test.alpha) | ...
                (tbl_mat_sorted.q1 < test.alpha) | ...
                (tbl_mat_sorted.q2 < test.alpha) | ...
                (tbl_mat_sorted.q3 < test.alpha);
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

input.tableCaption = 'Significant treatment differences';
input.tableLabel = 'te_direct';
input.tableBorders = 0;

input.data = tex.data;
input.tableRowLabels = tex.rows; % tex.fID_filt; % 
input.tableColLabels = tex.cols;

input.tableColumnAlignment = 'c';
input.dataFormat = {'%.3f',6};

latex = latexTable(input);


%% function for latex formatting

function uptr = upwards_trend(test)
    uptr = cell(test.nb_features,1);
    for i = 1:test.nb_features
        mean_oxt = mean(test.sample{i}(:,1),'omitnan');
        mean_plc = mean(test.sample{i}(:,2),'omitnan');
        if mean_oxt >= mean_plc
            sign = '$+$';
        else
            sign = '$-$';
        end
        uptr{i} = sign;
    end
end
