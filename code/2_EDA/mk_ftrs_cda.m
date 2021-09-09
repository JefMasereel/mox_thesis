close all
clear

%% CDA-based feature estimations for ElectroDermal Activity (EDA)

% using Ledalab package for EDA decomposition
% git clone https://github.com/ledalab/ledalab.git
% documentation and research papers on http://ledalab.de

% NOTE
% CDA is best for features that require high time resolution
% but limited in information on SCR morphology


%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\data\';
dir_org = [root 'ORIGINAL\rawdata_signals\'];
dir_eda = [root '2_EDA\'];
dir_cda = [root '2_EDA\leda_CDA\'];

% files to be included
filenames = load([dir_eda 'adjusted_filenames.mat']).filenames;


%% Some initial params

% fixed sample rate for EDA data
fs = 32; % [Hz]

% method params from Ledalab
batchmode = load([dir_cda 'batchmode_protocol.mat']);
methods.Leda_CDA = batchmode.protocol;
SCR_ampthreshold = methods.Leda_CDA.command.export_scrlist(1);


%% Specify features to be included in estimations

% directly derived features
ftr_names = {'Phasic_SC_Power','Diff2_SC_Power','SC_slope','Amp_sum'};

% (avg,std,mad) for SCR peak amplitude & inter-response-interval (IRI)
ftrs.aggreg = {'avg','std','mad'};
ftrs.scr_series = {'Amp','IRI'};   
for i = 1:length(ftrs.scr_series)
    for j = 1:length(ftrs.aggreg)
        ftr_names{end+1} = [ftrs.scr_series{i} '_' ftrs.aggreg{j}];
    end
end

% array for collection of results
ftr_values = zeros(length(filenames),length(ftr_names));


%% Collect feature estimations for each file

for i=1:length(filenames)
    
    disp(filenames{i})
    
    % load results from Ledalab
    CDA = load([dir_cda filenames{i} '.mat']);
    CDA_scr = load([dir_cda filenames{i} '_scrlist.mat']);
    
    % precalc polyfit for EDA slope
    EDA = CDA.data.conductance;
    coefs = polyfit(1:length(EDA),EDA,1);
    EDA_slope = coefs(1);
    
    % collect relevant results
    F1 = bandpower(CDA.analysis.phasicData);
    F2 = bandpower(diff(CDA.data.conductance,2));
    F3 = EDA_slope;
    F4 = CDA_scr.scrList.CDA.amp;
    F5 = diff(CDA_scr.scrList.CDA.onset);
    
    % remove NaN values from serial estimations for cleaner results
    % only if that leaves no samples left will ftr_values return NaN
    F4 = F4(~isnan(F4));
    F5 = F5(~isnan(F5));
    
    % collect aggregate features and assign to ftr_values matrix
    ftr_values(i,contains(ftr_names,'Phasic_SC_Power')) = F1;
    ftr_values(i,contains(ftr_names,'Diff2_SC_Power')) = F2;
    ftr_values(i,contains(ftr_names,'SC_slope')) = F3;
    ftr_values(i,contains(ftr_names,'Amp_sum')) = sum(F4);
    ftr_values(i,contains(ftr_names,'Amp_avg')) = mean(F4);
    ftr_values(i,contains(ftr_names,'Amp_std')) = std(F4);
    ftr_values(i,contains(ftr_names,'Amp_mad')) = mad(F4);
    ftr_values(i,contains(ftr_names,'IRI_avg')) = mean(F5);
    ftr_values(i,contains(ftr_names,'IRI_std')) = std(F5);
    ftr_values(i,contains(ftr_names,'IRI_mad')) = mad(F5);
end


%% add 'CDA_' prefix to all ftr names to distinguish from DDA ftrs

for i =1:length(ftr_names)
    ftr_names{i} = ['CDA_' ftr_names{i}];
end


%% Save logdata.mat

save([dir_eda 'logdata_cda.mat'],...
        'filenames',...
        'ftr_names',...
        'ftr_values',...
        'methods');


