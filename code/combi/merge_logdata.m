close all
clear

%% Feature analysis prep - merge logdata.mat files
% logdata_x files contain all feature estimations from each mode

% Script based on previous efforts in qssr_code repo
% as part of thesis by Jef Masereel


%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\data\';
dir_bvp = [root '1_BVP\'];
dir_eda = [root '2_EDA\'];
dir_rsp = [root '3_RSP\'];

% include all files where possible
% some will be missing for certain features
filenames = load([root 'adjusted_filenames.mat']).filenames;


%% collect logdata and merge into overarching logfile

log_prv = load([dir_bvp 'logdata_prv.mat']);
log_cda = load([dir_eda 'logdata_cda.mat']);
log_dda = load([dir_eda 'logdata_dda.mat']);
% log_rspa = load([dir_rsp 'alternative\logdata_rspa.mat']);
log_rspc = load([dir_rsp 'findpeaks\logdata_rspc.mat']);

% check if any logfiles need filename padding
assert(length(filenames)==length(log_cda.filenames)); % ok
assert(length(filenames)==length(log_dda.filenames)); % ok
assert(length(filenames)==length(log_prv.filenames)); % ok
% assert(length(filenames)==length(log_rspa.filenames));
assert(length(filenames)==length(log_rspc.filenames)); % ok

% pad PRV ftrs with nan where files are missing
padval_prv = nan(length(filenames),length(log_prv.ftr_names));
for i = 1:length(filenames)
    findmid = cellfun(@(c) isequal(c,filenames{i}), log_prv.filenames);
    if any(findmid)
        mid = find(findmid==1);
        padval_prv(i,:) = log_prv.ftr_values(mid,:);
    else
        % mid missing for this mode, leave nan value
    end
end

% add 'PRV_' prefix to PRV features
for i =1:length(log_prv.ftr_names)
    log_prv.ftr_names{i} = ['PRV_' log_prv.ftr_names{i}];
end

% replace 'RSPC_' prefix with 'RSP_' (when ignoring RSPA)
log_rspc.ftr_names = strrep(log_rspc.ftr_names,'RSPC_','RSP_');

% concatenate ftrs as [PRV CDA DDA RSP]
ftr_names = cat(2,log_prv.ftr_names',log_cda.ftr_names,log_dda.ftr_names,log_rspc.ftr_names); % removed log_rspa.ftr_names,
ftr_values = cat(2,padval_prv,log_cda.ftr_values,log_dda.ftr_values,log_rspc.ftr_values);     % removed log_rspa.ftr_values,


%% Collect metadata from method params

methods.Leda_CDA = log_cda.methods.Leda_CDA;
methods.Leda_DDA = log_dda.methods.Leda_DDA;
methods.Clas_PRV = log_prv.methods;
% methods.RSP_A = log_rspa.methods; % RSPA not included in final report
methods.RSP = log_rspc.methods;     % referred to as RSP in further reports


%% Save logdata_all.mat

save([root 'logdata_all.mat'],...
        'filenames',...
        'ftr_names',...
        'ftr_values',...
        'methods');