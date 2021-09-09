close all
clear

%% Feature analysis - build inspection plots for c52 (TE)
% 8 features, 3 conditions, compare oxt/plc directly as in te_compact.tex


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


%% Specify source & figure layout

% main test specs
test.paired = false;
test.method = 'ranksum'; %'ttest2';
test.alpha = .05;
verbose = false;

% features of interest (from te_compact.mat)
% ignore less sign. power ftr (HF) for even eight
ftrlist = {'DDA_HalfRec_time_avg',...
           'DDA_P50_width_avg'   ,...
           'DDA_P50_width_std'   ,...
           'DDA_P50_width_mad'   ,...
           'DDA_HalfRec_time_mad',...
           'DDA_Amp_sum'         ,...
           'DDA_HalfRec_time_std',...
           'PRV_HFn'             };
%            'PRV HF'              }; 
ftrx = get_idx(log,ftrlist);

ftrH = [true false true;  % ftr1
        true false true;  % ftr2
        true false false;
        true false false;
        true false false;
        true false true;
        true false true;
        false false true]; % ftr8

% fig specs
figrows = ftrlist;
figcols = {'rest','stress','dbcsr'};
nb_rows = length(figrows);
nb_cols = length(figcols);

% build sample structure
samples = cell(nb_rows,nb_cols);
sign_diff = ftrH;




%% Treatment effect - dftr,RS,OT/PL
% do the features find significant differences in treatment effects?

gethere = 'rest';
test.data = dftr.data;
test.idx = [(didx.RS & didx.oxytocin),...
            (didx.RS & didx.placebo)];
% build
test = build_samples(test,pids,dftr.rownames);
for r = 1:nb_rows                        % use nb_cols if transposed
    c = find(contains(figcols,gethere)); % use figrows if transposed
    samples{r,c} = test.sample{ftrx{r}}; % use c if transposed 
end


%% Treatment effect - dftr,EYE,OT/PL
% do the features find significant differences in treatment effects?

gethere = 'stress';
test.data   = dftr.data;
test.idx = [(didx.EYE & didx.oxytocin),...
            (didx.EYE & didx.placebo)];
% build
test = build_samples(test,pids,dftr.rownames);
for r = 1:nb_rows                        % use nb_cols if transposed
    c = find(contains(figcols,gethere)); % use figrows if transposed
    samples{r,c} = test.sample{ftrx{r}}; % use c if transposed 
end


%% Compare BCSR for treatment groups - DBCSR, oxt/plc

gethere = 'dbcsr';
test.data   = dbcsr.data;
test.idx = [dridx.oxytocin,...
            dridx.placebo];
% build
test = build_samples(test,pids,dbcsr.rownames);
for r = 1:nb_rows                        % use nb_cols if transposed
    c = find(contains(figcols,gethere)); % use figrows if transposed
    samples{r,c} = test.sample{ftrx{r}}; % use c if transposed 
end


%% build figures from results
% ugly.. move to python for violinplots? 

% %% violinplot
% figure
% for r = 1:nb_rows
%     for c = 1:nb_cols
%         subplot(nb_rows,nb_cols,(r-1)*nb_cols+c);
%         violinplot(samples{r,c});
%     end
% end
% 
% 
% %% boxplot & parallelcoords
% figure    
% coordLineStyle = 'k.';
% for r = 1:nb_rows
%     for c = 1:nb_cols
%         subplot(nb_rows,nb_cols,(r-1)*nb_cols+c);
%         boxplot(samples{r,c}, 'Symbol', coordLineStyle); hold on;
%         parallelcoords(samples{r,c}, 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
%           'Marker', '.', 'MarkerSize', 10);
%     end
% end
% 
% 
% %% only one feature
% figure    
% coordLineStyle = 'k.';
% for r = 1:nb_rows
% %     for c = 1:nb_cols
%         subplot(1,nb_rows,r);
%         boxplot(samples{r,1}, 'Symbol', coordLineStyle); hold on;
%         parallelcoords(samples{r,c}, 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
%           'Marker', '.', 'MarkerSize', 10);
% %     end
% end


%% get direction of trends
trends = cell(nb_rows,nb_cols);
for r = 1:nb_rows
    for c = 1:nb_cols
        mean_oxt = mean(samples{r,c}(:,1),'omitnan');
        mean_plc = mean(samples{r,c}(:,2),'omitnan');
        trends{r,c} = mean_oxt >= mean_plc;
    end
end

%% errorbar plots (compact)



figure
for r = 1:nb_rows
    for c = 1:nb_cols
        mean_oxt = mean(samples{r,c}(:,1),'omitnan');
        mean_plc = mean(samples{r,c}(:,2),'omitnan');
        std_oxt = std(samples{r,c}(:,1),'omitnan');
        std_plc = std(samples{r,c}(:,2),'omitnan');
        
        subplot(nb_rows,nb_cols,(r-1)*nb_cols+c);
        hold on;
        errorbar([mean_oxt,mean_plc],[std_oxt,std_plc],'x');
        if sign_diff(r,c) % add red line plot between means
            plot([1,2],[mean_oxt,mean_plc],'r');
        end
        xlim([0.5,2.5]);
        xlabel(strrep(ftrlist{r},'_','\_'));
        %xlabel(figcols{c});
    end
end

for i = 1:3
    subplot(nb_rows,nb_cols,i)
    title(figcols{i});
end


%% save for later use

save([dir_intermediates 'te_boxplots.mat'],'samples');



%% localfunction

function idx = get_idx(log,ftrlist)
    idx = cell(size(ftrlist));
    for i = 1:length(ftrlist)
        idx{i} = find(contains(log.ftr_names,ftrlist{i}));
        if length(idx{i}) > 1 
            disp('unexpected feature ID in ftrlist'); 
        end 
    end
end