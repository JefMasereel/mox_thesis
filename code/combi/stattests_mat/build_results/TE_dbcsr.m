close all
clear

%% Feature analysis - apply tests and build result reports
% using TFDATA from build_tfdata.m (transformations applied)

% treatment effects DBCSR
% oxt/plc for DBCSR -> confirm treatment effect VS accustomization


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

% build data handling structures & effect/response data (dftr,bc_eye)
% write new copy of log.mat to use prepped DATA iso. raw features
log = orglog; log.ftr_values = DATA;
hdl = build_handlers(log,placebo,attachment,filenames);

pids = hdl.pids;
idx = hdl.idx;
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
test.data   = dbcsr.data;
test.paired = false;
test.method = 'ranksum'; %'ttest2';
test.alpha = .05;

% toggle output
verbose = true;  % print assumption check results
plotFDR = true; % individual ranked pFDR plots
FDR_all = true;  % grouped figure with ranked results


%% Prepare table construction

% working with ftr values loaded from DATA
% apply 3 tests, get pvals & pFDR from each
% sort rows by pFDR ranking, only include if p<alpha

nb_features = size(test.data,2);

% table layout
tbl.varnames = {'p','q'};
tbl.rownames = log.ftr_names;


%% Compare BCSR for treatment groups - DBCSR, oxt/plc

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
tbl.pval = test.p;
tbl.qval = test.q;
tbl.pFDR = test.FDR;


%% Review results in matlab table

tbl_mat = table(tbl.pval,tbl.qval,...
                'VariableNames',tbl.varnames,...
                'RowNames',tbl.rownames);

% sort by first col: pval(pre)
[tbl_mat_sorted,I] = sortrows(tbl_mat); 

% only include significant rejections
filter_sorted = (tbl_mat_sorted.p < test.alpha) | ...
                (tbl_mat_sorted.q < test.alpha);
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

input.tableCaption = 'te_dbcsr';
input.tableLabel = 'te_dbcsr';
input.tableBorders = 0;

input.data = tex.data;
input.tableRowLabels = tex.rows; % ftr_ref.ID_str(filter);
input.tableColLabels = tex.cols;

input.tableColumnAlignment = 'c';
input.dataFormat = {'%.4f',2};

latex = latexTable(input);
