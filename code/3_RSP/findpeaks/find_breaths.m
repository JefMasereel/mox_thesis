close all
clear

%% RSP - Method C.1
% experimental design based on previous versions
% as part of thesis by Jef Masereel

% Step 1
% - basic preprocessing
% - annotate peaks (exhale onsets) with robust method
% - annotate corresponding troughs (inhale onsets)
% - translate TTP pairs to breath annotations
% - collect features for breath annotations
% - save to logdata_rsp_C.mat for analysis
% - save overview of annotations to .TIF

% ref RDECO by Jonathan Moeyersons:
% https://physionet.org/content/r-deco/1.0.0/
% git clone https://gitlab.esat.kuleuven.be/biomed-public/r-deco.git


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\data\';
dir_rsp = [root '3_RSP\findpeaks\'];
dir_org = [root 'ORIGINAL\rawdata_signals\'];

% list of included files (ref init_process.m)
% NOTE: some files were dropped during manual cleaning!
filenames = load([root 'filenames.mat']).filenames;


%% check for missing filenames and correct iterable

adj_files = {};
for i = 1:length(filenames)
    try 
        tmp = readtable([dir_org filenames{i} '.txt']);
        adj_files{end+1} = filenames{i};
    catch
        % if file is missing, filename is not added to adj_files list
    end
end
fprintf("%d missing files in dataset, removed from further analysis \n",...
    length(filenames)-length(adj_files));
disp(setdiff(filenames,adj_files))

% overwrite for now 
% possible to add other filling rules later if needed
filenames = adj_files;
save([dir_rsp 'adjusted_filenames.mat'],'filenames');


%% Preprocess RSP signals
% note: experimental method & parameters, cfr. thesis report

% general processing
fs_org = 128; % [Hz]
BP_window = [0.08 0.6]; 
fixlength = 300*fs_org;

% standardization method (optional, only relevant if RSP was calibrated)
StdzMethod = 'mean_mad';
apply_stdz = false;

% customized peak detection
QuantileRange = [0.30 0.70]; % use smaller range for data with many peaks
MinPromFactor = 0.40;        % use smaller factor to increase sensitivity

% breaths = zeros(size(filenames));
% features = zeros(size(filenames));
for i = 1:length(filenames)
    
    printable_filename = strrep(filenames{i},'_',' ');
    disp(printable_filename)
    
    tmp = readtable([dir_org filenames{i} '.txt']);
    rsp_raw = tmp.Var4';
    
    % fix length, rmv NaNs, BP filter
    rsp = basic_prep(rsp_raw,fs_org,fixlength,BP_window);

    % standardize amplitude for consistency in findpeaks args
    if apply_stdz
        rsp = custom_standardization(rsp,StdzMethod);
    end
    
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

    % save results as .tif file for later review
    savedir = [dir_rsp 'visualchecks\tif_annotations\'];
    print_results(breaths(i),filenames{i},savedir);
    
    % collect ftr results
    features(i) = ftr_breaths(breaths(i));

end


%% Reformat and save to logdata_rsp_C
% C suffix used to identify annotation method (cfr. 1_RSP/alternative/)

biomarkers = features(1).biomarkers;
scorenames = features(1).scorenames;

% collect feature values
nb_features = length(biomarkers)*length(scorenames);
ftr_values = zeros(length(filenames),nb_features);
ftr_names = cell(1,nb_features);
for i = 1:length(biomarkers)      % {'Ti','Te','RR','Vt','MV'}
    for j = 1:length(scorenames)  % {'MN','SD','MAD','CV','AR'}
        
        % define filling pattern
        col = (i-1)*length(scorenames) + j;
        
        % build variable names
        ftr_names{col} = ['RSPC_' biomarkers{i} '_' scorenames{j}];
        
        % collect data values
        for k = 1:length(filenames)
            ftr_values(k,col) = features(k).(scorenames{j}).(biomarkers{i});
        end
    end
end


% document applied methods
methods.prep = struct('fs_org',fs_org,'fs_ds',fs_org,...
                      'BP_window',BP_window,'fixlength',fixlength);
methods.standardization = struct('stdz_method',StdzMethod);
methods.findpeaks = struct('QuantileRange',QuantileRange,...
                           'MinPromFactor',MinPromFactor);


%% SAVE RESULTS

savepath_logdata = [dir_rsp 'logdata_rspc.mat'];

% Save & load logdata.mat
save(savepath_logdata,...
        'filenames',...
        'ftr_names',...
        'ftr_values',...
        'methods');
    
logdata = load(savepath_logdata);
                       

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


function print_results(breaths,filename,savedir)
% save .tif file with overview of annotation results

    inh = breaths.inh_onsets;           % location indices          [samples] 
    exh = breaths.exh_onsets;           % location indices          [samples] 
    rsp = breaths.respiration;          % standardized ampl @ fs    [a.u.]
    time = breaths.time;                % time signal @ fs          [sec]

    printable_filename = strrep(filename,'_',' ');
    figure('visible','off'); 
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 50 25], 'PaperUnits', 'centimeters', 'PaperSize', [50 25]);

    hold on
    plot(time,rsp,'b')
    plot(time(inh),rsp(inh),'g*')
    plot(time(exh),rsp(exh),'r*')
    legend('RSP, BP [0.08 0.6], standardized','inhale onsets','exhale onsets')
    title(['Annotation results method C, for file ' printable_filename]);
    xlabel('Time [sec]')
    ylabel('RSP amplitude [a.u.]')

    % save figure as .tif
    filepath = [savedir filename];
    print(gcf,filepath,'-dtiff','-r300');
end
