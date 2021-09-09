close all
clear

%% Feature analysis - collect features of interest
% treatment effects on activity in stress conditions


%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\';
datadir = [root 'data\'];

% save/load intermediate results in data folder
dir_intermediates = [datadir 'combi\intermediates_analysis\'];

% utility functions
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


%% General specs

% inspection: use as loaded from SOURCE
inspect.normality = SOURCE.normality;
inspect.vartestn = SOURCE.vartestn;

% rejection threshold
test.alpha = .05;

% toggle output
verbose = true;  % print assumption check results
plotFDR = false; % individual ranked pFDR plots
FDR_all = false;  % grouped figure with ranked results


%% Treatment effect - dftr,EYE,OT/PL
% do the features find significant differences in treatment effects?

% main test specs
test.paired = false;
test.method = 'ranksum'; %'ttest2';
test.data   = dftr.data;
test.idx    = [(didx.EYE & didx.oxytocin),...
               (didx.EYE & didx.placebo)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,dftr.rownames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features of interest:
TE_d = (test.h==1);


%% check each treatment group to confirm OT-specific effects


%% TE_ftrs for OT group

% test specs OT group
test.data   = DATA;
test.paired = true;
test.method = 'signrank'; %'ttest';
test.idx = [(idx.EYE & idx.oxytocin & idx.PRE),...
            (idx.EYE & idx.oxytocin & idx.POST)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features that capture treatment effects in the oxytocin group
TE_oxt = (test.h==1);
results.oxt = test.sample;


%% TE_ftrs for PL group

% test specs PL group
test.data   = DATA;
test.paired = true;
test.method = 'signrank'; %'ttest';
test.idx = [(idx.EYE & idx.placebo & idx.PRE),...
            (idx.EYE & idx.placebo & idx.POST)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features that capture treatment effects in the placebo group
TE_plc = (test.h==1);
results.plc = test.sample;


%% get features of interest by logical interpretation of results
% note: strict rejection based on alpha! can be more nuanced..

% "OT-specific effects"
OT_spec_effects = (TE_d==1 & TE_oxt==1 & TE_plc==0);
fprintf('\nOT-specific effect on EYE activity in %d features\n',...
    sum(OT_spec_effects));
for ftr = 1:length(log.ftr_names)
    if OT_spec_effects(ftr)
        fprintf('%s\n',...
            log.ftr_names{ftr});
    end
end

% "PL-specific effects"
PL_spec_effects = (TE_d==1 & TE_oxt==0 & TE_plc==1);
fprintf('\nPL-specific effect on EYE activity in %d features\n',...
    sum(PL_spec_effects));
for ftr = 1:length(log.ftr_names)
    if PL_spec_effects(ftr)
        fprintf('%s\n',...
            log.ftr_names{ftr});
    end
end

% shared effects (placebo, habituation, other)
nonspec_effects = (TE_d==0 & TE_oxt==1 & TE_plc==1);
fprintf('\nNonspecific effect on EYE activity in %d features\n',...
    sum(nonspec_effects));
for ftr = 1:length(log.ftr_names)
    if nonspec_effects(ftr)
        fprintf('%s\n',...
            log.ftr_names{ftr});
    end
end


%% check direction of OT-specific effects with barplots
% cfr /build_plots/te_boxplots.m

fprintf('\nDirection of OT-specific effects:\n');
for ftr = 1:length(log.ftr_names)
    if OT_spec_effects(ftr)
        mean_pre  = mean(results.oxt{ftr}(:,1));
        mean_post = mean(results.oxt{ftr}(:,2));
        if mean_pre > mean_post
            fprintf('pre- to post-treatment decrease in %s \n',log.ftr_names{ftr});
        else
            fprintf('pre- to post-treatment increase in %s \n',log.ftr_names{ftr});
        end
    end
end

fprintf('\nDirection of PL-specific effects:\n');
for ftr = 1:length(log.ftr_names)
    if PL_spec_effects(ftr)
        mean_pre  = mean(results.plc{ftr}(:,1));
        mean_post = mean(results.plc{ftr}(:,2));
        if mean_pre > mean_post
            fprintf('pre- to post-treatment decrease in %s \n',log.ftr_names{ftr});
        else
            fprintf('pre- to post-treatment increase in %s \n',log.ftr_names{ftr});
        end
    end
end

fprintf('\nDirection of nonspecific effects in OT group:\n');
for ftr = 1:length(log.ftr_names)
    if nonspec_effects(ftr)
        mean_pre  = mean(results.oxt{ftr}(:,1));
        mean_post = mean(results.oxt{ftr}(:,2));
        if mean_pre > mean_post
            fprintf('pre- to post-treatment decrease in %s \n',log.ftr_names{ftr});
        else
            fprintf('pre- to post-treatment increase in %s \n',log.ftr_names{ftr});
        end
    end
end

fprintf('\nDirection of nonspecific effects in PL group:\n');
for ftr = 1:length(log.ftr_names)
    if nonspec_effects(ftr)
        mean_pre  = mean(results.plc{ftr}(:,1));
        mean_post = mean(results.plc{ftr}(:,2));
        if mean_pre > mean_post
            fprintf('pre- to post-treatment decrease in %s \n',log.ftr_names{ftr});
        else
            fprintf('pre- to post-treatment increase in %s \n',log.ftr_names{ftr});
        end
    end
end
