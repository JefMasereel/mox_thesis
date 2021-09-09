close all
clear

%% BVP - Classical PRV analysis, feature estimation

% Script based on previous efforts in qssr_code repo
% as part of thesis by Jef Masereel

% interpretation of fiducial points based on prior work by J. Lazaro
% ref J. Lazaro paper: DOI 10.1109/JBHI.2013.2267096


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\';

% save logdata_prv at dir_bvp\
dir_hdr = [root 'data\'];
dir_bvp = [root 'data\1_BVP\'];
dir_original  = [dir_hdr 'ORIGINAL\rawdata_signals\'];

% choose nM annotations as fiducial points
dir_validfids = [dir_hdr '1_BVP\valid_fids\'];
dir_fid = [dir_validfids 'vnM\'];
fid_suffix = '_vnM.mat';

% path to save .tif files for later review
dir_tif = [dir_bvp 'tif_prv\'];

% list of included files (ref mk_validfids.m)
% NOTE: some files were dropped during manual cleaning! 
tmp = load([dir_validfids 'adjusted_filenames.mat']);
filenames = tmp.filenames;
clear tmp

% use all files now (possible to remove filenames if wanted)
ignore_msrs = [];
included_msrs = setdiff(filenames, ignore_msrs);


%% Define method params

% standardized segment length [seconds]
Lsec = 300;

% reconstruct PRV in desired format
fs_prv = 4;         %Hz, target fs for spectral analysis
fs_fid = 128;       %Hz, sample rate during fiducial calculation
type = 'flattened'; % flattened ends improve PSD estimation

% Apply BP filter to focus on relevant spectral content of PRV
N  = 3;                    % filter order
Wn = [0.03 1.00]*2/fs_prv; % norm. cutoff frequency
[B,A] = butter(N,Wn,'bandpass');

% compute Welch PSD of PRV
windowsize = 40; % [s]
overlap = 50;    % [%]
dfft_pts = 2^10; % [-]

% border frequencies for classical PRV
freq_Ld = 0.04; % [Hz]
freq_Lu = 0.15; % [Hz]
freq_Hd = 0.15; % [Hz]
freq_Hu = 0.40; % [Hz]
freq_TOd = 0.00; % [Hz]
freq_TOu = 1.00; % [Hz]
border_freqs = [freq_Ld, freq_Lu; ...   % LF range
                freq_Hd, freq_Hu; ...   % HF range
                freq_TOd, freq_TOu];    % range for total power

% document applied methods
ftr_names = {'PR_avg';'PR_std';'PR_mad';'TOT';'LF';'HF';'LFn';'HFn';'LH_ratio'};

methods.make_prv = struct(...
    'Lsec',Lsec,'fs_prv',fs_prv,'fs_fid',fs_fid,'type',type);

methods.filtfilt = struct(...
    'N',N,'Wn',Wn,'B',B,'A',A);

methods.pwelch = struct(...
    'windowsize',windowsize,'overlap',overlap,'dfft_pts',dfft_pts);

methods.get_PRV_classic = struct(...
    'border_freqs',border_freqs);


%% Define data structure to handle inputs and outputs

% count size of iterables
nb_valid_msrs = size(filenames,1) - size(ignore_msrs,1);
nb_valid_ftrs = size(ftr_names,1);

% data matrix contents
ftr_values = zeros(nb_valid_msrs, nb_valid_ftrs);


%% Process data

for m = 1:length(included_msrs)
    
    disp(included_msrs{m})
    printable_filename = strrep(included_msrs{m},'_',' ');
    
    % get fid data
    tmp = load([dir_fid included_msrs{m} fid_suffix]);
    fiducials = tmp.vnM;
    
    % make (flattened!) PRV
    PRV = make_prv(Lsec,fs_prv,fs_fid,fiducials,type);

    % apply BP filter
    PRV_filt = filtfilt(B,A,PRV);

    % compute Welch PSD
    [PSD,f] = pwelch(PRV_filt,...
                   fs_prv*windowsize,...
                   fs_prv*windowsize*overlap/100,...
                   dfft_pts,...
                   fs_prv);

    % compute classical power contents
    results = get_PRV_classic(PSD,f,border_freqs);
    
    % collect direct BVP features (PR)
    pulse_rate = diff(fiducials)/fs_fid;
    ftr_values(m,1) = mean(pulse_rate);
    ftr_values(m,2) = std(pulse_rate);
    ftr_values(m,3) = mad(pulse_rate);
    
    % collect classical PRV features
    for f = 4:nb_valid_ftrs
        ftr_values(m,f) = results.(ftr_names{f});
    end
    
    % save figures for later review
%     figure('visible','off'); 
%     set(gcf, 'Units', 'centimeters', 'Position', [0 0 12 5], 'PaperUnits', 'centimeters', 'PaperSize', [12 5]);
%     time = 1/fs_prv:1/fs_prv:length(PRV)/fs_prv;
%     plot(time,PRV)
%     title(['PRV for ' printable_filename])
%     xlabel('time [seconds]')
%     ylabel('PRV [a.u.]')
%     filepath = [dir_tif included_msrs{m}];
%     print(gcf,filepath,'-dtiff','-r300');
end


%% Save data structure to local dir

% for consistency in further processing
filenames = included_msrs;

save([dir_bvp 'logdata_prv.mat'],...
        'filenames',...
        'ftr_names',...
        'ftr_values',...
        'methods');
    

%% Local functions

function PRV = make_prv(Lsec,fs_prv, fs_fid,fiducials,type)
% This function builds the PRV signal to specified fs
% from RDECO savefile format
% to interpolated time signal at fs_prv sample freq

% WARNING : this script assumes that all time signal have been cut to a
% fixed length in previous processing steps (before deriving fiducials).
% Make_prv will return invalid results if Lsec does not match the data. 

% Lsec          [time]      fixed length of input signal, in seconds
% fs_prv        [Hz]        target fs for pulse rate variability signal
% fs_fid        [Hz]        fs used for calculation of fiducial points
% fiducials     [samples]   fiducial points from earlier annotations

% RR_int [time]         duration between fid points, in samples @fs_fid
% R_loc  [time]         location of fid points, in matlab time format

% type          'flat' for flattened ends (better PSD results)
%               'normal' for direct PRV, non-flattened
% PRV_flat      PRV signal @fs_prv, with flattened ends

% use spline to build PRV from fiducial points
x  = fiducials(2:end);                          % [samples]
y  = diff(fiducials)/fs_fid*1000;               % [ms]
xq = 0:fs_fid/fs_prv:((Lsec*fs_fid)-1);
PRV = spline(x,y,xq);
t_prv = 0:1/fs_prv:((Lsec*fs_fid)-1)/fs_fid;

% flatten the signal ends to prevent artefacts in spectral analysis
if type=='flattened'
    PRV_flat = PRV;
    % flatten all known PRV values before first measurement
    i1 = find(t_prv<fiducials(2)/fs_fid,1,'last');
    PRV_flat(1:i1) = PRV(i1);
    % flatten all known PRV values after last measurement
    in = find(t_prv>fiducials(end)/fs_fid,1,'first');
    PRV_flat(in:end) = PRV(in);
    % return flattened PRV
    PRV = PRV_flat;
end
end


function results = get_PRV_classic(PSD,f,border_freqs)
% PSD, f       : spectral power content of given time signal
% border_freqs : define low and high frequency ranges (LF & HF) 
%  [Ld, Lu ;   : lower and upper range for LF band
%   Hd, Hu ;   : lower and upper range for HF band
%   TOd, TOu]  : lower and upper range for total band of interest

% results struct with fields:
%   LF, HF      absolute power in defined spectral ranges
%   TOT         absolute total power in given spectral data
%   LFn, HFn    normalized power in defined spectral ranges

% sanity check on border_freqs size
assert(all(size(border_freqs)==[3,2]),'invalid border_freqs param, wrong size')

% find nearest discrete border indices
didx = zeros(3,2);
for i = 1:3
    didx(i,1) = find(f<=border_freqs(i,1),1,'last');  % lower bound
    didx(i,2) = find(f>=border_freqs(i,2),1,'first'); % upper bound
end

% sanity check, indices should be within available bounds
assert(~any(isempty(didx)),'invalid border_freqs param, idx out of bounds')

% compute power within defined freq ranges
LF  = trapz(f(didx(1,1):didx(1,2)),PSD(didx(1,1):didx(1,2)));
HF  = trapz(f(didx(2,1):didx(2,2)),PSD(didx(2,1):didx(2,2)));
TOT = trapz(f(didx(3,1):didx(3,2)),PSD(didx(3,1):didx(3,2)));

% normalized values
LFn = LF/TOT;
HFn = HF/TOT;

% LF-HF ratio
LH_ratio = LFn/HFn;

% gather in struct for easier reference
results = struct('TOT',TOT,...             % total power
                 'LF',LF,...               % abs power in LF
                 'HF',HF,...               % abs power in HF
                 'LFn',LFn,...             % rel power in LF
                 'HFn',HFn,...             % rel power in HF
                 'LH_ratio',LH_ratio);     % power ratio
end

