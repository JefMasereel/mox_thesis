close all
clear

%% saam scatter plots, based on saam_te & saam_pre


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

saam.scores = attachment.scores(:,3:5);
saam.varnames = attachment.names(3:5);

% build data handling structures & treatment effect data (dftr)
% write new copy of log.mat to use prepped DATA iso. raw features
log = orglog; log.ftr_values = DATA;
hdl = build_handlers(log,placebo,attachment,filenames);
pids = hdl.pids;
idx = hdl.idx;
aidx = hdl.aidx;
ridx = hdl.ridx;
bcsr = hdl.bcsr;
dftr = hdl.dftr;
didx = hdl.didx;
dridx = hdl.dridx;
dbcsr = hdl.dbcsr;

% general specs
test.alpha = .05;
test.paired = true;
test.method = 'spearman_corr';

        
%% RS PRE - figure 1

%           SAAM scale  SAAM scale      feature ID
pairlist = {'rest'      'security',     'RSP_RR_MN';
            'stress',   'anxiety',      'CDA_IRI_std'};
nb_plots = size(pairlist,1);
test.data = DATA;

figure();

% subplot 1
plot_id = 1;
test.idx  = (idx.PRE & idx.RS);
saam.var = pairlist{plot_id,2};
test = build_samples(test,pids,filenames);
test = attach_saam(test,saam,pids);
test = apply_maintest(test);
ftrx = find(contains(log.ftr_names,pairlist{plot_id,3}));

subplot(1,nb_plots,plot_id)
scatter(test.sample{ftrx}(:,1),test.sample{ftrx}(:,2),'filled');
xlabel([strrep(pairlist{plot_id,3},'_',' ') '  (' pairlist{plot_id,1} ')']);
ylabel(pairlist{plot_id,2});
title(['\rho = ' num2str(test.rho(ftrx),3)]);

% subplot 2
plot_id = 2;
test.idx  = (idx.PRE & idx.EYE);
saam.var = pairlist{plot_id,2};
test = build_samples(test,pids,filenames);
test = attach_saam(test,saam,pids);
test = apply_maintest(test);
ftrx = find(contains(log.ftr_names,pairlist{plot_id,3}));

subplot(1,nb_plots,plot_id)
scatter(test.sample{ftrx}(:,1),test.sample{ftrx}(:,2),'filled');
xlabel([strrep(pairlist{plot_id,3},'_',' ') '  (' pairlist{plot_id,1} ')']);
ylabel(pairlist{plot_id,2});
title(['\rho = ' num2str(test.rho(ftrx),3)]);



%% SAAM TE - figure 2


%           TE      condition   SAAM scale  feature ID
pairlist = {'oxt'    'stress'    'security', 'DDA_SC_slope';
            'plc',   'rest',     'anxiety',  'DDA_HalfRec_time_std';
            'plc',   'stress',   'anxiety',  'CDA_Amp_avg';
            'plc',   'stress',   'security', 'PRV_LF'};
nb_plots = size(pairlist,1);

test.data = dftr.data;
figure();

% subplot 1
plot_id = 1;
test.idx  = (didx.EYE & didx.oxytocin);

saam.var = pairlist{plot_id,3};
test = build_samples(test,pids,dftr.rownames);
test = attach_saam(test,saam,pids);
test = apply_maintest(test);
ftrx = find(contains(log.ftr_names,pairlist{plot_id,4}));
subplot(1,nb_plots,plot_id)
scatter(test.sample{ftrx}(:,1),test.sample{ftrx}(:,2),'filled');
xlabel([strrep(pairlist{plot_id,4},'_',' ') ' (' pairlist{plot_id,1} ', ' pairlist{plot_id,2} ')']);
ylabel(pairlist{plot_id,3});
title(['\rho = ' num2str(test.rho(ftrx),3)]);


% subplot 2
plot_id = 2;
test.idx  = (didx.RS & didx.placebo);

saam.var = pairlist{plot_id,3};
test = build_samples(test,pids,dftr.rownames);
test = attach_saam(test,saam,pids);
test = apply_maintest(test);
ftrx = find(contains(log.ftr_names,pairlist{plot_id,4}));
subplot(1,nb_plots,plot_id)
scatter(test.sample{ftrx}(:,1),test.sample{ftrx}(:,2),'filled');
xlabel([strrep(pairlist{plot_id,4},'_',' ') ' (' pairlist{plot_id,1} ', ' pairlist{plot_id,2} ')']);
ylabel(pairlist{plot_id,3});
title(['\rho = ' num2str(test.rho(ftrx),3)]);


% subplot 3
plot_id = 3;
test.idx  = (didx.EYE & didx.placebo);

saam.var = pairlist{plot_id,3};
test = build_samples(test,pids,dftr.rownames);
test = attach_saam(test,saam,pids);
test = apply_maintest(test);
ftrx = find(contains(log.ftr_names,pairlist{plot_id,4}));
subplot(1,nb_plots,plot_id)
scatter(test.sample{ftrx}(:,1),test.sample{ftrx}(:,2),'filled');
xlabel([strrep(pairlist{plot_id,4},'_',' ') ' (' pairlist{plot_id,1} ', ' pairlist{plot_id,2} ')']);
ylabel(pairlist{plot_id,3});
title(['\rho = ' num2str(test.rho(ftrx),3)]);


% subplot 4
plot_id = 4;
test.idx  = (didx.EYE & didx.placebo);

saam.var = pairlist{plot_id,3};
test = build_samples(test,pids,dftr.rownames);
test = attach_saam(test,saam,pids);
test = apply_maintest(test);
ftrx = find(contains(log.ftr_names,pairlist{plot_id,4}),1,'first');
subplot(1,nb_plots,plot_id)
scatter(test.sample{ftrx}(:,1),test.sample{ftrx}(:,2),'filled');
xlabel([strrep(pairlist{plot_id,4},'_',' ') ' (' pairlist{plot_id,1} ', ' pairlist{plot_id,2} ')']);
ylabel(pairlist{plot_id,3});
title(['\rho = ' num2str(test.rho(ftrx),3)]);



