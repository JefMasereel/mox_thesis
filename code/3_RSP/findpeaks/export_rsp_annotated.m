close all
clear


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\data\';
dir_rsp = [root '3_RSP\findpeaks\'];
dir_from = [root 'ORIGINAL\rawdata_signals\'];
dir_to = [dir_rsp 'export_rsp\annotated\'];

% list of included files (ref init_process.m)
% NOTE: some files were dropped during manual cleaning!
tmp = load([dir_rsp 'adjusted_filenames.mat']);
filenames = tmp.filenames;
clear tmp


%% Preprocess RSP signals

% general processing
fs_org = 128; % [Hz]
BP_window = [0.08 0.6];
fixlength = 300*fs_org;

% standardization method
StdzMethod = 'mean_mad';

% customized peak detection
QuantileRange = [0.30 0.70]; % use smaller range for data with many peaks
MinPromFactor = 0.40;        % use smaller factor to increase sensitivity

% breaths = zeros(size(filenames));
% features = zeros(size(filenames));
for i = 1:length(filenames)
    
    printable_filename = strrep(filenames{i},'_',' ');
    disp(printable_filename)
    
    tmp = readtable([dir_from filenames{i} '.txt']);
    rsp_raw = tmp.Var4';
    
    % fix length, rmv NaNs, BP filter
    rsp = basic_prep(rsp_raw,fs_org,fixlength,BP_window);

    % preliminary ampl stdz disabled to avoid adding uncertain variations
%     % standardize amplitude for consistency in findpeaks args
%     rsp = custom_standardization(rsp,StdzMethod);
    
    % peak detection threshold
    amp_swing = diff(quantile(rsp,QuantileRange));
    MinPeakProminence = MinPromFactor*amp_swing;
    
    % apply findpeaks with robust args
    [peak_amp,peak_idx] = findpeaks(rsp,'MinPeakProminence',MinPeakProminence);
    
    % collect troughs corresponding to each peak
    trough_amp = zeros(1,length(peak_amp)-1);
    trough_idx = zeros(1,length(peak_idx)-1);
    for k=2:length(peak_idx)
        [amp,idx] = min(rsp(peak_idx(k-1):peak_idx(k)));
        trough_amp(k-1) = amp;
        trough_idx(k-1) = peak_idx(k-1)+idx-1;
    end
    % ignore first exhale onset peak for consistent annotation sequence
    % makes it a lot easier to compute breath metrics for paired TTP
    % TTP pairs always inhale to exhale onset (trough to peak)
    peak_amp(1) = []; peak_idx(1) = [];
    assert(all(size(peak_idx)==size(trough_idx)));
    
    % interpretation of annotations:
    breaths(i).respiration = rsp;
    breaths(i).srate = fs_org;
    breaths(i).inh_onsets = trough_idx;
    breaths(i).exh_onsets = peak_idx;
    breaths(i).time = 1/fs_org:1/fs_org:length(rsp)/fs_org;

    % export annotated respiration to .mat for external use
    respiration = breaths(i);
    savepath = [dir_to filenames{i} '.mat'];
    save(savepath,'respiration');
end


%% Local functions

function rsp = basic_prep(rsp,fs,fixlength,BP_window)
    
    % Set fixed length @ 300 seconds
    try
        rsp = rsp(1:fixlength);
    catch
        fprintf('signal has only %d samples \n',length(rsp));
    end

    % remove nan values at end
    rsp = rsp(~isnan(rsp)); 
    
    % bandpass filter
    N  = 3; 
    Wn = BP_window*2/fs;
    [B,A] = butter(N,Wn,'bandpass');
    rsp = filtfilt(B,A,rsp);   
end

function std_sig = custom_standardization(sig,method)
% adjust standardization method to tune sensitivity to large breaths
    
    if strcmp(method,'zscore')
        std_sig = (sig-mean(sig))/std(sig);
    elseif strcmp(method,'med_mad')
        std_sig = (sig-median(sig))/mad(sig);
    elseif strcmp(method,'mean_mad')
        std_sig = (sig-mean(sig))/mad(sig);
    else
        fprintf("Warning - unknown method specified for signal standardization\n")
    end
end