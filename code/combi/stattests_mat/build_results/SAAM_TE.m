close all
clear

%% Feature analysis - apply tests and build result reports
% using TFDATA from build_tfdata.m (transformations applied)

% attachment score interactions (C_spearman) (pvals,pFDR)
% with treatment effects for each group

% for i={anxiety, avoidance, security}
% 1. corr(SAAM(i),ftr(TE,rs,oxt)
% 2. corr(SAAM(i),ftr(TE,rs,plc)
% 3. corr(SAAM(i),ftr(TE,eye,oxt)
% 4. corr(SAAM(i),ftr(TE,eye,plc)
% 5. corr(SAAM(i),ftr(TE,eye-rs,oxt)
% 6. corr(SAAM(i),ftr(TE,eye-rs,plc)


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
aidx = hdl.aidx;
didx = hdl.didx;
dftr = hdl.dftr;
dridx = hdl.dridx;
dbcsr = hdl.dbcsr;


%% Notes on test and sample specification in this script
% applied some OOP-ish variable construction to clean up testing process

% for correlations, use additional SAAM info
% > test    = build_samples(..)           build first sample as usual
% > test    = attach_saam(saam,var)       attach 1 SAAM score as 2nd sample
% > inspect = inspect_distribution(test)  continue as usual
% > test    = apply_maintest

% Basic construction as with difference tests
% test.paired = true;
% test.method = 'ttest';
% test.data   = log.ftr_values;
% test.idx  = (idx.PRE & idx.RS);    % sample A from ftr data

% Additional SAAM score info
% saam.scores = attachment.scores(:,3:5);
% saam.varnames = attachment.names(3:5);
% sample B from SAAM score data is matched to sample A by subject ID

% note: idx arrays should contain logical values!
% using numerical ones and zeros breaks the indexing..

% apply_maintest(method)
% pearson_corr    parametric      paired    correlation
% spearman_corr   non-parametric  paired    correlation


%% General specs

% inspection: use as loaded from SOURCE
inspect.normality = SOURCE.normality;
inspect.vartestn = SOURCE.vartestn;

% main test specs
test.paired = true;
test.method = 'spearman_corr';
test.alpha = .05;

% correlation tests require additional attachment information
saam.scores = attachment.scores(:,3:5);
saam.varnames = attachment.names(3:5);

% toggle output
verbose = true;  % print assumption check results
plotFDR = false; % individual ranked pFDR plots
FDR_all = true;  % grouped figure with ranked results


%% Prepare collection of results

% 3x3x2 = (test,saam,treatment)
tbl.pval = cell(3,3,2);
tbl.qval = cell(3,3,2);
tbl.pFDR = cell(3,3,2);
tbl.rho = cell(3,3,2);
tbl.rownames = log.ftr_names;


%% SAAM correlation - dftr,rs,oxt
% do the SAAM scores correlate with inter-subject variations
% after treatment with oxytocin?

test.data = dftr.data;
test.idx  = (didx.RS & didx.oxytocin);

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,dftr.rownames);
    test    = attach_saam(test,saam,pids);
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
    tbl.pval{1,i,1} = test.p;
    tbl.qval{1,i,1} = test.q;
    tbl.pFDR{1,i,1} = test.FDR;
    tbl.rho{1,i,1}  = test.rho;
end


%% SAAM correlation - dftr,rs,plc
% do the SAAM scores correlate with inter-subject variations
% after treatment with placebo?

test.data = dftr.data;
test.idx  = (didx.RS & didx.placebo);

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,dftr.rownames);
    test    = attach_saam(test,saam,pids);
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
    tbl.pval{1,i,2} = test.p;
    tbl.qval{1,i,2} = test.q;
    tbl.pFDR{1,i,2} = test.FDR;
    tbl.rho{1,i,2}  = test.rho;
end


%% SAAM correlation - dftr,eye,oxt
% do the SAAM scores correlate with inter-subject variations
% after treatment with oxytocin?

test.data = dftr.data;
test.idx  = (didx.EYE & didx.oxytocin);

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,dftr.rownames);
    test    = attach_saam(test,saam,pids);
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
    tbl.pval{2,i,1} = test.p;
    tbl.qval{2,i,1} = test.q;
    tbl.pFDR{2,i,1} = test.FDR;
    tbl.rho{2,i,1}  = test.rho;
end


%% SAAM correlation - dftr,eye,plc
% do the SAAM scores correlate with inter-subject variations
% after treatment with placebo?

test.data = dftr.data;
test.idx  = (didx.EYE & didx.placebo);

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,dftr.rownames);
    test    = attach_saam(test,saam,pids);
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
    tbl.pval{2,i,2} = test.p;
    tbl.qval{2,i,2} = test.q;
    tbl.pFDR{2,i,2} = test.FDR;
    tbl.rho{2,i,2}  = test.rho;
end


%% SAAM correlation - DBCSR,oxt

test.data = dbcsr.data;
test.idx  = dridx.oxytocin;

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,dbcsr.rownames);
    test    = attach_saam(test,saam,pids);
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
    tbl.pval{3,i,1} = test.p;
    tbl.qval{3,i,1} = test.q;
    tbl.pFDR{3,i,1} = test.FDR;
    tbl.rho{3,i,1}  = test.rho;
end


%% SAAM correlation - DBCSR,plc

test.data = dbcsr.data;
test.idx  = dridx.placebo;

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,dbcsr.rownames);
    test    = attach_saam(test,saam,pids);
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
    tbl.pval{3,i,2} = test.p;
    tbl.qval{3,i,2} = test.q;
    tbl.pFDR{3,i,2} = test.FDR;
    tbl.rho{3,i,2}  = test.rho;
end


%% Review results in matlab tables (one for each treatment group)

%% saam_te, oxytocin
% compressed format compared to other scripts due to size of results

tbl.rownames = strrep(pad(log.ftr_names),'_',' ');
tbl.varnames.headers = {'feature','sample','SAAM','pval','qval'};
tbl.varnames.test = {'rest','stress','bc-stress'};
tbl.varnames.saam = saam.varnames; % {'anx','sec','avd'};

testnames = {'te_oxytocin','te_placebo'};
for t=1:2           % treatment group
    
    fprintf('\n%s\n',testnames{t});
    
    for i = 1:3         % sample / test number

    %     fprintf('\n');
        fprintf('\\hline \n');

        for j = 1:3     % SAAM vars (anx,sec,avd)

            [tmp.pval,I] = sort(tbl.pval{i,j,t});
            filter = find(tmp.pval < test.alpha); %% replace with qval ??

%             tmp.pFDR = tbl.pFDR{i,j,t}(I);
            tmp.qval = tbl.qval{i,j,t}(I);
            tmp.rho  = tbl.rho{i,j,t}(I);
            tmp.ftr  = tbl.rownames(I);
            tmp.test = tbl.varnames.test{i};
            tmp.saam = tbl.varnames.saam{j};

            for k = 1:length(filter)

                fprintf('%s & %s & %s & %.4f & %.4f & %.4f %s \n',...
                    tmp.ftr{k},tmp.test,tmp.saam,...
                    tmp.pval(k),tmp.qval(k),tmp.rho(k),'\\');
            end
        end
    end
end


%% OUTDATED below here

% %% Table 1, p-values
% 
% tbl_mat = table(tbl.pval{1,1},tbl.pval{2,1},tbl.pval{3,1},tbl.pval{4,1},...
%                 tbl.pval{1,2},tbl.pval{2,2},tbl.pval{3,2},tbl.pval{4,2},...
%                 tbl.pval{1,3},tbl.pval{2,3},tbl.pval{3,3},tbl.pval{4,3},...
%                 'VariableNames',tbl.varnames(1,:),...
%                 'RowNames',tbl.rownames);
% 
% % sort by first col: pval(pre)
% [tbl_mat_sorted,I] = sortrows(tbl_mat); 
% 
% % only include significant rejections
% filter_sorted = (tbl_mat_sorted.p1_RS < test.alpha) | ...
%                 (tbl_mat_sorted.p2_RS < test.alpha) | ...
%                 (tbl_mat_sorted.p3_RS < test.alpha) | ...
%                 (tbl_mat_sorted.p4_RS < test.alpha) | ...
%                 (tbl_mat_sorted.p1_EYE < test.alpha) | ...
%                 (tbl_mat_sorted.p2_EYE < test.alpha) | ...
%                 (tbl_mat_sorted.p3_EYE < test.alpha) | ...
%                 (tbl_mat_sorted.p4_EYE < test.alpha);
% tbl_tex = tbl_mat_sorted(filter_sorted,:);
% 
% if verbose
%     disp(tbl_tex)
% end
% 
% % ftrIDs sorted & filtered for compact view
% tex.fID_sort = ftr_ref.ID_str(I);
% tex.fID_filt = ftr_ref.ID_str(filter_sorted);
% 
% % printable header names
% tex.rows = strrep(tbl_tex.Properties.RowNames,'_',' ');
% tex.cols = tbl_tex.Properties.VariableNames;
% tex.data = tbl_tex.Variables;
% 
% 
% % build latex table
% fprintf('Latex table\n\n');
% 
% input.data = tex.data;
% input.tableRowLabels = tex.rows; % ftr_ref.ID_str(filter);
% input.tableColLabels = tex.cols;
% 
% input.tableBorders = 0;
% input.tableColumnAlignment = 'r';
% input.dataFormat = {'%.4f',6};
% 
% latex = latexTable(input);


%% Visual inspection SAAM score distributions 
% draft here, might be good to make proper figure in python

figure
plotmatrix(saam.scores(:,:))
title('SAAM score distributions')
