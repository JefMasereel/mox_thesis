close all
clear

%% Feature analysis - collect features of interest
% treatment effects on baseline corrected stress response (bcsr)


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


%% Treatment effect - dbcsr,OT/PL
% do the features find significant differences in treatment effects?

% main test specs
test.paired = false;
test.method = 'ranksum'; %'ttest2';
test.data   = dbcsr.data;
test.idx    = [dridx.oxytocin,...
               dridx.placebo];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,dbcsr.rownames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features of interest:
TE_d = (test.h==1);


%% check each treatment group to confirm OT-specific effects


%% TE_bcsr for OT group

% test specs OT group
test.paired = true;
test.method = 'signrank';
test.data   = bcsr.data;
test.idx    = [(ridx.PRE  & ridx.oxytocin),...
               (ridx.POST & ridx.oxytocin)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,bcsr.rownames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features that capture treatment effects in the oxytocin group
TE_oxt = (test.h==1);
results.oxt = test.sample;


%% TE_bcsr for PL group

% test specs PL group
test.paired = true;
test.method = 'signrank';
test.data   = bcsr.data;
test.idx    = [(ridx.PRE  & ridx.placebo),...
               (ridx.POST & ridx.placebo)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,bcsr.rownames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features that capture treatment effects in the oxytocin group
TE_plc = (test.h==1);
results.plc = test.sample;


%% get features of interest by logical interpretation of results
% note: strict rejection based on alpha! can be more nuanced..

% "OT-specific effects"
OT_spec_effects = (TE_d==1 & TE_oxt==1 & TE_plc==0);
fprintf('\nOT-specific effect on BCSR in %d features\n',...
    sum(OT_spec_effects));
for ftr = 1:length(log.ftr_names)
    if OT_spec_effects(ftr)
        fprintf('%s\n',...
            log.ftr_names{ftr});
    end
end

% "PL-specific effects"
PL_spec_effects = (TE_d==1 & TE_oxt==0 & TE_plc==1);
fprintf('\nPL-specific effect on BCSR in %d features\n',...
    sum(PL_spec_effects));
for ftr = 1:length(log.ftr_names)
    if PL_spec_effects(ftr)
        fprintf('%s\n',...
            log.ftr_names{ftr});
    end
end

% shared effects (placebo, habituation, other)
nonspec_effects = (TE_d==0 & TE_oxt==1 & TE_plc==1);
fprintf('\nNonspecific effect on BCSR in %d features\n',...
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
            fprintf('pre- to post-treatment decrease in SR amplitude of %s \n',log.ftr_names{ftr});
        else
            fprintf('pre- to post-treatment increase in SR amplitude of %s \n',log.ftr_names{ftr});
        end
    end
end

fprintf('\nDirection of PL-specific effects:\n');
for ftr = 1:length(log.ftr_names)
    if PL_spec_effects(ftr)
        mean_pre  = mean(results.plc{ftr}(:,1));
        mean_post = mean(results.plc{ftr}(:,2));
        if mean_pre > mean_post
            fprintf('pre- to post-treatment decrease in SR amplitude of %s \n',log.ftr_names{ftr});
        else
            fprintf('pre- to post-treatment increase in SR amplitude of %s \n',log.ftr_names{ftr});
        end
    end
end

fprintf('\nDirection of nonspecific effects in OT group:\n');
for ftr = 1:length(log.ftr_names)
    if nonspec_effects(ftr)
        mean_pre  = mean(results.oxt{ftr}(:,1));
        mean_post = mean(results.oxt{ftr}(:,2));
        if mean_pre > mean_post
            fprintf('pre- to post-treatment decrease in SR amplitude of %s \n',log.ftr_names{ftr});
        else
            fprintf('pre- to post-treatment increase in SR amplitude of %s \n',log.ftr_names{ftr});
        end
    end
end

fprintf('\nDirection of nonspecific effects in PL group:\n');
for ftr = 1:length(log.ftr_names)
    if nonspec_effects(ftr)
        mean_pre  = mean(results.plc{ftr}(:,1));
        mean_post = mean(results.plc{ftr}(:,2));
        if mean_pre > mean_post
            fprintf('pre- to post-treatment decrease in SR amplitude of %s \n',log.ftr_names{ftr});
        else
            fprintf('pre- to post-treatment increase in SR amplitude of %s \n',log.ftr_names{ftr});
        end
    end
end


%% check origin of stress response (cfr. diff_sr.m)
% the previous tests indicate the significant reductions in SR
% the following takes a closer look at the significance of SRs pre/post


%% RS-EYE difference OT, PRE

% main test specs
test.data   = DATA;
test.paired = true;
test.method = 'signrank';
test.idx    = [(idx.PRE & idx.oxytocin & idx.RS),...
               (idx.PRE & idx.oxytocin & idx.EYE)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features that capture a response in the OT group, PRE treatment
SR_oxt_pre = (test.h==1);
results.oxt_pre = test.sample;


%% RS-EYE difference OT, POST

% main test specs
test.data   = DATA;
test.paired = true;
test.method = 'signrank';
test.idx    = [(idx.POST & idx.oxytocin & idx.RS),...
               (idx.POST & idx.oxytocin & idx.EYE)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features that capture a response in the OT group, POST treatment
SR_oxt_post = (test.h==1);
results.oxt_post = test.sample;


%% RS-EYE difference PL, PRE

% main test specs
test.data   = DATA;
test.paired = true;
test.method = 'signrank';
test.idx    = [(idx.PRE & idx.placebo & idx.RS),...
               (idx.PRE & idx.placebo & idx.EYE)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features that capture a response in the PL group, PRE treatment
SR_plc_pre = (test.h==1);
results.plc_pre = test.sample;


%% RS-EYE difference PL, POST

% main test specs
test.data   = DATA;
test.paired = true;
test.method = 'signrank';
test.idx    = [(idx.POST & idx.placebo & idx.RS),...
               (idx.POST & idx.placebo & idx.EYE)];

% check validity of assumptions, then apply test
test    = build_samples(test,pids,filenames);
inspect = inspect_distribution(test,inspect,verbose);
test    = apply_maintest(test);
[test.FDR,test.q] = mafdr(test.p);

% features that capture a response in the PL group, POST treatment
SR_plc_post = (test.h==1);
results.plc_post = test.sample;


%% get features of interest by logical interpretation of results
% note: strict rejection based on alpha! can be more nuanced..

% OT-specific changes in amplitude/significance of stress response
% OT_SR_toggle = (TE_oxt==1 & SR_oxt_pre~=SR_oxt_post & TE_d==1 & TE_plc==0);
OT_SR_toggle = (TE_oxt==1 & SR_oxt_pre~=SR_oxt_post);
fprintf('\nOT-specific BCSR toggle in %d features\n',sum(OT_SR_toggle));

% PL-specific changes in amplitude/significance of stress response
PL_SR_toggle = (TE_plc==1 & SR_plc_pre~=SR_plc_post);
fprintf('PL-specific BCSR toggle in %d features\n',sum(PL_SR_toggle));

% check direction of toggles (binary switch)
fprintf('\nDirection of OT-specific BCSR pre/post toggles:\n');
for ftr = 1:length(log.ftr_names)
    if OT_SR_toggle(ftr)
        if SR_oxt_pre(ftr)==1 && SR_oxt_post(ftr)==0
            fprintf('SR toggles from h=1 to h=0 for %s \n',log.ftr_names{ftr});
        else
            fprintf('SR toggles from h=0 to h=1 for %s \n',log.ftr_names{ftr});
        end
    end
end
