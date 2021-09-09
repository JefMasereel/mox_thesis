close all
clear

%% Feature analysis - apply tests and build result reports
% using TFDATA from build_tfdata.m (transformations applied)

% attachment score interactions (C_spearman) (pvals,pFDR)
% with baseline activity at rest & in stress response

% for i={anxiety, avoidance, security}
% 1. corr(SAAM(i),ftr(pre,rs,all)
% 2. corr(SAAM(i),ftr(pre,eye,all)
% 3. corr(SAAM(i),ftr(pre,eye-rs,all)


%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\';
datadir = [root 'data\'];

% save/load intermediate results in data folderthyujjj
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
ridx = hdl.ridx;
bcsr = hdl.bcsr;


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


%% Prep datastructure for collection of results

% 3 tests x 3 SAAM scores
tbl.pval = cell(3,3);
tbl.qval = cell(3,3);
tbl.pFDR = cell(3,3);
tbl.rho  = cell(3,3);


%% SAAM correlation - pre,rs,all
% do the SAAM scores correlate with inter-subject variations before treatment?

test.data = DATA;
test.idx  = (idx.PRE & idx.RS);

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,filenames);
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
    tbl.pval{1,i} = test.p;
    tbl.qval{1,i} = test.q;
    tbl.pFDR{1,i} = test.FDR;
    tbl.rho{1,i} = test.rho;
end


%% SAAM correlation - pre,eye,all
% do the SAAM scores correlate with inter-subject variations before treatment?

test.data = DATA;
test.idx  = (idx.PRE & idx.EYE);

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,filenames);
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
    tbl.pval{2,i} = test.p;
    tbl.qval{2,i} = test.q;
    tbl.pFDR{2,i} = test.FDR;
    tbl.rho{2,i} = test.rho;
end


%% SAAM correlation - pre,bcsr,all
% do the SAAM scores correlate with inter-subject variations before treatment?

test.data = bcsr.data;
test.idx  = ridx.PRE;

for i = 1:length(saam.varnames)
    saam.var = saam.varnames{i};

    % check validity of assumptions, then apply test
    test    = build_samples(test,pids,bcsr.rownames);
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
    tbl.pval{3,i} = test.p;
    tbl.qval{3,i} = test.q;
    tbl.pFDR{3,i} = test.FDR;
    tbl.rho{3,i} = test.rho;
end


%% Review results in summary table 
% compressed format compared to other scripts due to size of results

tbl.rownames = strrep(pad(log.ftr_names),'_',' ');
tbl.varnames.headers = {'feature','sample','SAAM','p','q','rho'};
tbl.varnames.test = {'rest','stress','bc-stress'};
tbl.varnames.saam = saam.varnames; % {'anx','sec','avd'};

for i = 1:3         % sample / test number
    
%     fprintf('\n');
        fprintf('\\hline \n');
    
    for j = 1:3     % SAAM vars (anx,sec,avd)
        
        [tmp.pval,I] = sort(tbl.pval{i,j});
        filter = find(tmp.pval < test.alpha);  %% replace with qval ??
        
%         tmp.pFDR = tbl.pFDR{i,j}(I);
        tmp.qval = tbl.qval{i,j}(I);
        tmp.rho  = tbl.rho{i,j}(I);
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


%% Visual inspection SAAM score distributions 
% draft here, might be good to make proper figure in python

figure
[H,AX,BigAx,P,PAx] = plotmatrix(saam.scores(:,:));

scatter_range = [0.5,7.5];

for x=1:3
    for y=1:3
        AX(x,y).XLim = scatter_range;
        AX(x,y).YLim = scatter_range;
    end
end

for x=1:3
    xlabel(AX(3,x),saam.varnames{x})
end

for y=1:3
    ylabel(AX(y,1),saam.varnames{y})
end
