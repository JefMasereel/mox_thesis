close all
clear

%% BVP - collect cleaned nA annotations and derive full fiducials

% Script based on previous efforts in qssr_code repo
% as part of thesis by Jef Masereel

% interpretation of fiducial points based on prior work by J. Lazaro
% ref J. Lazaro paper: DOI 10.1109/JBHI.2013.2267096


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\';

dir_hdr = [root 'data\'];
dir_bvp = [root 'data\1_BVP\'];
dir_original  = [dir_hdr 'ORIGINAL\rawdata_signals\'];
dir_rdeco_out = [dir_bvp 'rdeco_out\'];

dir_validfids = [dir_hdr '1_BVP\valid_fids\'];
dir_vnA = [dir_validfids 'vnA\'];
dir_vnB = [dir_validfids 'vnB\'];
dir_vnM = [dir_validfids 'vnM\'];

% list of included files
tmp = load([dir_bvp 'adjusted_filenames.mat']);
filenames = tmp.filenames;
clear tmp

% add path to some custom local functions
addpath([root 'code\1_BVP\localfunctions']);


%% Process parameters

% Start: assume rnA is already cleaned in RDECO
% rnA --> vnA : set chosen fs, tune annotations based on local max
% vnA --> vnB : set window [ms]
% vnB --> vnM : automatic

% downsample signal from fs_org to fs_int
% or keep original fs (fs_int=nan)
fs_org = 128;       % [Hz]
fs_int = nan;       % [Hz]

% set standardized length for all measurements
fixed_time = 300;   % [s] Default 300s

% set window for nA tuning
window_nA = 0.020;  % [s] Default 20ms (keep small!)

% set search window for nB estimation
window_nB = 0.150;  % [s] Default 150ms (make larger if necessary)


%% Some basic precalcs

if isnan(fs_int)
    fs = fs_org;
else
    fs = fs_int;
end

fixed_size  = fixed_time*fs;
windowsize_nA = get_valid_windowsize(window_nA,fs,true);
windowsize_nB = get_valid_windowsize(window_nB,fs,true);


%% check for missing filenames and correct iterable

adj_files = {};
for i = 1:length(filenames)
    try 
        tmp = load([dir_rdeco_out filenames{i} '_nA.mat']);
        adj_files{end+1} = filenames{i};
    catch
        % if file is missing, filename is not added to adj_files list
    end
end
fprintf("%d missing files in rdeco_out, removed from further analysis \n",...
    length(filenames)-length(adj_files));

% overwrite for now 
% possible to add other filling rules later if needed
filenames = adj_files;
save([dir_validfids 'adjusted_filenames.mat'],'filenames');


%% APPLY ANNOTATION METHODS - initialize vars

% % Remarks
% rnA data is stored as .mat files with data.[R_loc,RR_int] from RDECO
% R_loc is expressed in [seconds], use this to map onto chosen fs
% here it is assumed that rnA has already been cleaned in RDECO

% to review impact of tuning manual annotations
% Sum of Squared Errors (SSE) for nA
SSE_xA = zeros(size(filenames));


%% APPLY ANNOTATION METHODS - main loop
for i = 1:length(filenames)
    
%     disp(filenames{i})
    
    %% Get nA (assume cleaned already)

    % read corresponding time signal (PPG)
    sig = get_ppg(dir_original,filenames{i},fs_org,fs_int,fixed_size);
    
    % read rnA
    tmp = load([dir_rdeco_out filenames{i} '_nA.mat']);
    nA_idx = seconds(tmp.data.R_loc{1});                % [sec]
    nA_prv = tmp.data.RR_int{1};                        % [ms]
    
    % Convert annotations to samples (from seconds) @fs
    rnA = round(nA_idx*fs);
    % sanity checks (cfr ex51)
    check_rnA_formatting(sig,fs,rnA,nA_prv,fixed_size,filenames{i});
    
    % finetune manual annotations
    [vnA,SSE_xA(i)] = tune_rnA(rnA,sig,windowsize_nA,fixed_size);
    

    %% Get nB
    
    rnB = find_rnB(sig,fs,vnA,windowsize_nB);
    vnB = rnB;
    
    
    %% Get nM
    
    rnM = find_rnM(sig,vnA,vnB);
    vnM = rnM;
    

    %% Collect results
    
    % save results
    save([dir_vnA filenames{i} '_vnA.mat'],'vnA');
    save([dir_vnB filenames{i} '_vnB.mat'],'vnB');
    save([dir_vnM filenames{i} '_vnM.mat'],'vnM');

    % update progress file
%     progress(i,'clean')=1;
    
end