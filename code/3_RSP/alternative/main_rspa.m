close all
clear

%% RSP - Method A, automated artefact correction
% note: comparative discussion of methods can be found in thesis report

% uses BreathMetrics tool for breath annotations
% git clone https://github.com/zelanolab/breathmetrics.git
addpath('C:\Users\jefma\repos\breathmetrics\');
addpath('C:\Users\jefma\repos\breathmetrics\breathmetrics_functions\')


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\';
dir_hdr = [root 'data\'];
dir_rsp = [dir_hdr '3_RSP\alternative\'];
dir_org = [dir_hdr 'ORIGINAL\rawdata_signals\'];
dir_tif = [dir_rsp 'visualchecks\'];

% list of included files (ref init_process.m)
% NOTE: some files were dropped during manual cleaning! 
filenames = load([dir_hdr 'filenames.mat']).filenames;

% custom functions for artefact detection and masking
addpath([root 'code\3_RSP\alternative\localfunctions\']);


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


%% prep method variables

% general processing
fs_org = 128; % [Hz]
BP_window = [0.08 0.6];
fixlength = 300*fs_org;

% for BM process
zScore = 1;
verbose = 0;
simplify = 1;
baselineCorrectionMethod = 'sliding';

% for artefact correction
method.name = 'stdev';
method.rmoutliers = 'grubbs';
windowsize = 714;               % [samples] based on findings ex85
minsize_segment = 3*windowsize; % [samples] practical minimum length


%% apply artefact correction with full preprocessing

for i = 1:length(filenames)
    
    printable_filename = strrep(filenames{i},'_',' ');
    disp(printable_filename)
    
    tmp = readtable([dir_org filenames{i} '.txt']);
    rsp_raw = tmp.Var4';
    
    rsp = basic_prep(rsp_raw,fs_org,fixlength,BP_window);
    
    bmObj = breathmetrics(rsp, fs_org, 'humanBB');
    bmObj.correctRespirationToBaseline(baselineCorrectionMethod,zScore,verbose);
    
    mask = artf_mask(bmObj.baselineCorrectedRespiration,method,windowsize);
    segments_all = get_clean_segments(bmObj.baselineCorrectedRespiration,mask);
    segments = rmv_segments(segments_all,minsize_segment);
    
    % save .tif figures for review
    savedir = [dir_tif 'maskresults\'];
    print_maskresults(bmObj.baselineCorrectedRespiration,...
                      mask,fs_org,filenames{i},savedir);

    % reapply BM to each segment & collect results
    % do NOT reapply zscore standardization! baselinecorrection is ok.
    results = cell(size(segments.start));
    breaths = cell(size(segments.start));
    for k = 1:length(segments.start)
        rsp_k = segments.full{k};
        bmObj_k = breathmetrics(rsp_k,fs_org,'humanBB');
        bmObj_k.estimateAllFeatures(0,'sliding',1,0);
        results{k} = ftr_rsp(bmObj_k);
        breaths{k} = {bmObj_k.inhaleOnsets;bmObj_k.exhaleOnsets};
    end
    
    % merge feature scores across segments
    % merge breath annotations across segments as well
    mid_results = merge_results(results);
    merge_breaths(bmObj,breaths,segments);
    
    % save .tif figures incl. annotations for review
    savedir = [dir_tif 'annotations\'];
    print_BMresults(bmObj,mask,filenames{i},savedir);
    
    % how much of the signal is left for analysis?
    ratio_removedsignal = rmv_ratio(rsp,segments);
    
    % collect results into collective structure
    features(i) = mid_results;
    features_add(i) = ratio_removedsignal;
end


%% Reformat and save to logdata_rsp

biomarkers = features(1).biomarkers;
scorenames = features(1).scorenames;

% save mid_results and add the rmvsignal ratios
nb_features = length(biomarkers)*length(scorenames) + 1;
ftr_values = zeros(length(filenames),nb_features);
ftr_names = cell(1,nb_features);

% collect mid_results
for i = 1:length(biomarkers)      % {'Ti','Te','RR','Vt','MV'}
    for j = 1:length(scorenames)  % {'MN','SD','MAD','CV','AR'}
        
        % define filling pattern
        col = (i-1)*length(scorenames) + j;
        
        % build variable names
        ftr_names{col} = ['RSPA_' biomarkers{i} '_' scorenames{j}];
        
        % collect data values
        for k = 1:length(filenames)
            ftr_values(k,col) = features(k).(scorenames{j}).(biomarkers{i});
        end
    end
end

% collect ratio data as additional feature
ftr_names{end} = 'rsp_used_full';
ftr_values(:,end) = features_add';

% document applied methods
methods.prep = struct('fs_org',fs_org,'fs_ds',fs_org,...
                      'BP_window',BP_window,'fixlength',fixlength);
methods.artf_corr = struct('metric',method.name,...
                           'windowsize',windowsize,...
                           'minsize',minsize_segment);
methods.rmoutliers = struct('args','none, all default');


%% SAVE RESULTS

savepath_logdata = [dir_rsp 'logdata_rspa.mat'];

% Save & load logdata.mat
save(savepath_logdata,...
        'filenames',...
        'ftr_names',...
        'ftr_values',...
        'methods');
logdata = load(savepath_logdata);


%% local functions

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


function print_maskresults(rsp,mask,fs,filename,savedir)
% local function to save figures as tif files

    printable_filename = strrep(filename,'_',' ');
    time = 1/fs:1/fs:length(rsp)/fs;
    figure('visible','off'); 

    % define figure dimensions [cm]
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 50 25], 'PaperUnits', 'centimeters', 'PaperSize', [50 25]);

    % build figure
    hold on
    plot(time,rsp)
    plot(time,mask)
    hold off
    xlabel('time [s]')
    ylabel('RSP amplitude [a.u.]')
    legend('processed RSP signal','masking signal')
    title(['artefact detection for ' printable_filename])

    % save figure as .tif
    filepath = [savedir filename];
    print(gcf,filepath,'-dtiff','-r300');
end


function segments = get_clean_segments(rsp,mask)

    % basic reference signal from first derivative
    diff_mask = diff(mask);

    % set correct indicator at first sample 
    if mask(1)==0
        diff_mask(1) = -1;                      % set start flag
    elseif mask(1)==1
        diff_mask(1) = 0;                       % avoid stop flag
    end

    % set correct indicator at last sample 
    if mask(end)==0
        diff_mask(end) = +1;                    % set stop flag
    elseif mask(end)==1
        diff_mask(1) = 0;                       % avoid start flag
    end

    % collect reference idx for start & stop of each clean segment 
    % (vice versa for artf segments, but not needed here)
    segments.start = find(diff_mask==-1);
    segments.stop  = find(diff_mask==+1);

    % lengths should match if edge cases were solved correctly
    try
        assert(length(segments.start)==length(segments.stop));
        assert(all(segments.start<segments.stop));
    catch
       fprintf('artefact correction - unexpected start/stop indices found \n')
    end
    
    % collect segments as separate rsp signals
    segments.full = cell(size(segments.start));
    for i = 1:length(segments.full)
        segments.full{i} = rsp(segments.start(i):segments.stop(i));
    end
end


function segments_fix = rmv_segments(segments,minsize)
    
    segments_fix.start = double.empty();
    segments_fix.stop  = double.empty();
    segments_fix.full  = cell.empty();

    % ignore segments shorter than minsize [samples]    
    for i = 1:length(segments.full)
        if segments.stop(i) - segments.start(i) + 1 > minsize
            segments_fix.start(end+1) = segments.start(i);
            segments_fix.stop(end+1)  = segments.stop(i);
            segments_fix.full{end+1}  = segments.full{i};
        end
    end
end


function ratio = rmv_ratio(rsp,segments)
    
    L_org = length(rsp);
    L_use = 0;
    for k = 1:length(segments.start)
        L_use = L_use + length(segments.full{k});
    end
    ratio = L_use/L_org;
end

function mid_results = merge_results(results)
% results in format specified by ftr_rsp.m

    ftr.scorenames = results{1}.scorenames;
    ftr.biomarkers = results{1}.biomarkers;
    for b=1:length(ftr.biomarkers)
        ftr.full.(ftr.biomarkers{b}) = double.empty();
    end
    
    % concatenate biomarker arrays from segments
    for k = 1:length(results)
        for b=1:length(ftr.biomarkers)
            aa = ftr.full.(ftr.biomarkers{b});
            bb = results{k}.full.(ftr.biomarkers{b});
            ftr.full.(ftr.biomarkers{b}) = [aa bb];
        end
    end
    
    % recompute scores for merged biomarkers
    for b=1:length(ftr.biomarkers)
        tmp = ftr.full.(ftr.biomarkers{b});

        if length(tmp) <= 1
            [MN,SD,MAD,CV,AR] = deal(NaN);
            warning = true;
            note = 'insufficient breaths detected in segment';
            fprintf('Warning - %d breath(s) detected in segment, skipping ftrs \n',length(tmp));
        else
            MN = mean(tmp);
            SD = std(tmp);
            MAD = mad(tmp);
            CV = SD/MN;
            AR = corr(tmp(1:end-1)',tmp(2:end)');
            warning = false;
            note = '';
        end

        ftr.MN.(ftr.biomarkers{b})  = MN;
        ftr.SD.(ftr.biomarkers{b})  = SD;
        ftr.MAD.(ftr.biomarkers{b}) = MAD;
        ftr.CV.(ftr.biomarkers{b})  = CV;
        ftr.AR.(ftr.biomarkers{b})  = AR;
        ftr.warning.(ftr.biomarkers{b}) = {warning;note};
    end
    
    % return merged results data
    mid_results = ftr;
end


function merge_breaths(bmObj,breaths,segments)
% recombine breath annotations for visual inspection
% assign resulting annotations to bmObj of mid
    
    inhaleOnsets = [];
    exhaleOnsets = [];

    for k = 1:length(segments.start)
        % idx in segment reference
        inhaleOnsets_k = breaths{k}{1};
        exhaleOnsets_k = breaths{k}{2};
        
        % translate idx to mid reference
        ref_idx = segments.start(k);
        inhaleOnsets = [inhaleOnsets inhaleOnsets_k+ref_idx];
        exhaleOnsets = [exhaleOnsets exhaleOnsets_k+ref_idx];
    end
    
    % assign to bmObj of mid
    bmObj.inhaleOnsets = inhaleOnsets;
    bmObj.exhaleOnsets = exhaleOnsets;
end


function print_BMresults(bmObj,mask,filename,savedir)
% save .tif files of full annotation results after masking
    
    % get relevant info from mid bmObj
    fs = bmObj.srate;
    time = bmObj.time;
    resp = bmObj.baselineCorrectedRespiration;
    inh_idx = bmObj.inhaleOnsets;
    exh_idx = bmObj.exhaleOnsets;

    printable_filename = strrep(filename,'_',' ');
    figure('visible','off'); 
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 50 25], 'PaperUnits', 'centimeters', 'PaperSize', [50 25]);

    % build figure
    hold on
    plot(time,resp)
    plot(time,mask)
    plot(time(inh_idx),resp(inh_idx),'r*')
    plot(time(exh_idx),resp(exh_idx),'b*')
    legend('Processed RSP','Masking flag','inhale onsets','exhale onsets')
    title(['Annotation results with artefact correction for ' printable_filename])
    
    % save figure as .tif
    filepath = [savedir filename];
    print(gcf,filepath,'-dtiff','-r300');
end



