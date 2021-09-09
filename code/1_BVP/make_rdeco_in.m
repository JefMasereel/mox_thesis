close all
clear

%% BVP - preparing datafiles for cleaning in RDECO

% ref RDECO by Jonathan Moeyersons:
% https://physionet.org/content/r-deco/1.0.0/


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\';

dir_hdr = [root 'data\'];
dir_bvp = [root 'data\1_BVP\'];
dir_original = [dir_hdr 'ORIGINAL\rawdata_signals\'];
dir_rdeco_in = [dir_bvp 'rdeco_in\'];

% path to save .tif files for later review
dir_tif = [dir_bvp 'tif_bvp\'];


%% load filenames to be included in analysis 
% build with init_process.m before starting signal analysis!

tmp = load([dir_hdr 'filenames.mat']);
filenames = tmp.filenames;
clear tmp


%% define preprocessing parameters

% downsample (optional)
% currently opted to maintain highest resolution available
fs_org = 128;
fs_dns = 128; 

% bandpass filter (optional)
% for this dataset, BVP signals already clean enough
BP_order = 3;
BP_range = [0.3 40];
[B,A] = butter(BP_order,BP_range*2/fs_org,'bandpass');
apply_BP = true;

% set fixed length (optional)
fixed_length = 300; % [s]
apply_fix = true;


%% check for missing filenames and correct iterable

adj_files = {};
for i = 1:length(filenames)
    try 
        tmp = readtable([dir_original filenames{i} '.txt']);
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
save([dir_bvp 'adjusted_filenames.mat'],'filenames');


%% iterate over full dataset and apply preprocessing

for i = 1:length(filenames)
    
    printable_filename = strrep(filenames{i},'_',' ');
    
    % BPV signal at third column of source data
    tmp = readtable([dir_original filenames{i} '.txt']);
    sig = tmp.Var3;
    
    %% remove nan values at end if present
    % 1 second tolerance at end (usually nans at last few samples)
    
    if any(find(isnan(sig))<length(sig)-1*fs_org) 
        fprintf('Warning - NaN values removed in signal body \n')
    end
    sig = sig(~isnan(sig));
    
    
    %% apply BP filter
    if apply_BP
        sig = filtfilt(B,A,sig);
    end

    %% downsample
    t = 0:1/fs_org:(length(sig)-1)/fs_org;  % time @fs_org
    ppg.t = t(1) : 1/fs_dns : t(end);       % time @fs_dns
    ppg.sig_preproc = spline(t, sig, ppg.t);
    ppg.fs = fs_dns;

    %% cut to fixed length
    if apply_fix
        ppg = fix_length(ppg,fixed_length);
    end

    %% control printouts
    disp(filenames{i})                    
    disp(length(sig)/fs_org)
    disp(length(ppg.sig_preproc)/ppg.fs)
    
    %% save ppg data under filename
    save([dir_rdeco_in filenames{i} '.mat'],'ppg')

    %% TIFF figures for signal inspection
%     figure('visible','on'); 
%     set(gcf, 'Units', 'centimeters', 'Position', [0 0 12 5], 'PaperUnits', 'centimeters', 'PaperSize', [12 5]);
%     plot(ppg.t,ppg.sig_preproc)
%     title(['BVP for ' printable_filename])
%     xlabel('time [seconds]')
%     ylabel('BVP [a.u.]')
%     filepath = [dir_tif filenames{i}];
%     print(gcf,filepath,'-dtiff','-r300');
end



%% Local function

function fix = fix_length(ppg,max_length)
% shorten all time series in given ppg struct to specified length [seconds]
% ppg.fs specifies the used sampling rate

    max_idx = max_length * ppg.fs; % [samples]
    
    if length(ppg.t) <= max_idx
        fix = ppg;
    else
        fix.fs = ppg.fs;
        fix.t = ppg.t(1:max_idx);
        fix.sig_preproc = ppg.sig_preproc(1:max_idx);
    end
end