close all
clear

%% DDA-based feature estimations for ElectroDermal Activity (EDA)
% latest version

% Script based on previous efforts in qssr_code repo
% as part of thesis by Jef Masereel

% using Ledalab package for EDA decomposition
% git clone https://github.com/ledalab/ledalab.git
% documentation and research papers on http://ledalab.de

% NOTE
% DDA provides needed info to quantify morphology of the SCRs
% but at cost of lower time resolution compared to DDA


%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\data\';
dir_org = [root 'ORIGINAL\rawdata_signals\'];
dir_eda = [root '2_EDA\'];
dir_dda = [root '2_EDA\leda_DDA\'];

% files to be included
filenames = load([dir_eda 'adjusted_filenames.mat']).filenames;


%% Some initial params

% fixed sample rate for EDA data
fs = 32; % [Hz]

% method params from Ledalab
batchmode = load([dir_dda 'batchmode_protocol.mat']);
methods.Leda_DDA = batchmode.protocol;
SCR_ampthreshold = methods.Leda_DDA.command.export_scrlist(1);


%% Specify features to be included in estimations

% directly derived features
ftr_names = {'Phasic_SC_Power','Diff2_SC_Power','SC_slope','Amp_sum','Area_sum'};

% collect (avg,std,mad) for available estimates
ftrs.aggreg = {'avg','std','mad'};
ftrs.scr_series = {'Amp','IRI','Area','HalfRec_time','P50_width','R50_slope','P50_slope'};   
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
    DDA = load([dir_dda filenames{i} '.mat']);
    DDA_scr = load([dir_dda filenames{i} '_scrlist.mat']);
    
    % precalc polyfit for EDA slope
    EDA = DDA.data.conductance;
    coefs = polyfit(1:length(EDA),EDA,1);
    EDA_slope = coefs(1);
    
    % compute shape features for DDA impulses
    form = shape_ftrs(DDA,DDA_scr);
    
    % collect relevant results
    F1 = bandpower(DDA.analysis.phasicData);
    F2 = bandpower(diff(DDA.data.conductance,2));
    F3 = EDA_slope;
    F4 = DDA_scr.scrList.DDA.amp;
    F5 = diff(DDA_scr.scrList.DDA.onset);
    F6 = form.SCR_areas;
    F7 = form.R50_time;
    F8 = form.P50_width;
    F9 = form.R50_slope;
    F10 = form.P50_slope;
    
    % remove NaN values from serial estimations for cleaner results
    % only if that leaves no samples left will ftr_values return NaN
    F4 = F4(~isnan(F4));    F5 = F5(~isnan(F5));
    F7 = F7(~isnan(F7));    F8 = F8(~isnan(F8));
%     F9 = F9(~isnan(F9));    F10 = F10(~isnan(F10));
    
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
    ftr_values(i,contains(ftr_names,'Area_avg')) = mean(F6);
    ftr_values(i,contains(ftr_names,'Area_std')) = std(F6);
    ftr_values(i,contains(ftr_names,'Area_mad')) = mad(F6);    
    ftr_values(i,contains(ftr_names,'Area_sum')) = sum(F6);    
    ftr_values(i,contains(ftr_names,'HalfRec_time_avg')) = mean(F7); % change to R50t?
    ftr_values(i,contains(ftr_names,'HalfRec_time_std')) = std(F7);
    ftr_values(i,contains(ftr_names,'HalfRec_time_mad')) = mad(F7);
    ftr_values(i,contains(ftr_names,'P50_width_avg')) = mean(F8);
    ftr_values(i,contains(ftr_names,'P50_width_std')) = std(F8);
    ftr_values(i,contains(ftr_names,'P50_width_mad')) = mad(F8);
    ftr_values(i,contains(ftr_names,'R50_slope_avg')) = mean(F9);
    ftr_values(i,contains(ftr_names,'R50_slope_std')) = std(F9);
    ftr_values(i,contains(ftr_names,'R50_slope_mad')) = mad(F9);
    ftr_values(i,contains(ftr_names,'P50_slope_avg')) = mean(F10);
    ftr_values(i,contains(ftr_names,'P50_slope_std')) = std(F10);
    ftr_values(i,contains(ftr_names,'P50_slope_mad')) = mad(F10);
end


%% add 'DDA_' prefix to all ftr names to distinguish from CDA ftrs

for i =1:length(ftr_names)
    ftr_names{i} = ['DDA_' ftr_names{i}];
end


%% Save logdata.mat

save([dir_eda 'logdata_dda.mat'],...
        'filenames',...
        'ftr_names',...
        'ftr_values',...
        'methods');


%% Local function

function form = shape_ftrs(DDA,DDA_scr)
% compute series of half-recovery times, width percentiles and slopes
% for all impulse estimations corresponding to SCRs above threshold 0.05

% find registered SCRs (above amp_threshold)
ref_onset = DDA_scr.scrList.DDA.onset;
ref_idx = zeros(size(ref_onset));
for i = 1:length(ref_onset)
    ref_idx(i) = find(DDA.analysis.impulsePeakTime==ref_onset(i),1);
end

% % sanity check
% for i=1:length(ref_idx)
%     index = ref_idx(i);
%     disp(DDA_scr.scrList.DDA.amp(i)-DDA.analysis.amp(index))
% end

%% recovery features

% all estimations in [samples]
R50_time = zeros(size(ref_idx));    % halfrecovery time
R50_slope = zeros(size(ref_idx));   % interpolated slope of halfrecovery
P50_width = zeros(size(ref_idx));   % pulsewidth @ 50% of peak amplitude
P50_slope = zeros(size(ref_idx));   % direct slope from peak to P50
for i = 1:length(ref_idx)
    
    impulse = DDA.analysis.impulse{i};
    [peak_amp,peak_idx] = max(impulse);
    
    try
        R50_time(i) = find(impulse(peak_idx:end)<peak_amp*0.50,1);
        P50_width(i) = sum(impulse>peak_amp*0.50);
    catch
        fprintf('Warning - Impulse does not reach halfmax with 12s for impulse %d) \n',i);
        R50_time(i) = nan;
        P50_width(i) = nan;
    end
    
    % slopes are more robust to recoveries that don't reach 50% in 12s
    % just interpolate the available signal if recovery is too slow
    stop_idx = peak_idx + find(impulse(peak_idx:end)>peak_amp*0.50,1,'last') -1;

    R50_slope(i) = (impulse(stop_idx)-peak_amp)/stop_idx;
    coefs = polyfit(peak_idx:stop_idx,impulse(peak_idx:stop_idx),1);
    P50_slope(i) = coefs(1);
end

form.R50_time = R50_time;
form.R50_slope = R50_slope;
form.P50_width = P50_width;
form.P50_slope = P50_slope;

%% area features (estimated as area under curve of phasic signal)
form.DDA_areas = DDA.analysis.area;             % for all responses
form.SCR_areas = DDA.analysis.area(ref_idx);    % for SCR>amp_thresh

end


