close all
clear

%% BVP - review fiducial point annotations from 1_BVP\valid_fids

% Optional script for visual inspection of annotation results. 

% interpretation of fiducial points based on prior work by J. Lazaro
% ref J. Lazaro paper: DOI 10.1109/JBHI.2013.2267096


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\';

dir_hdr = [root 'data\'];
dir_bvp = [root 'data\1_BVP\'];
dir_original  = [dir_hdr 'ORIGINAL\rawdata_signals\'];
dir_rdeco_out = [dir_bvp '\rdeco_out\'];

dir_validfids = [dir_hdr '1_BVP\valid_fids\'];
dir_vnA = [dir_validfids 'vnA\'];
dir_vnB = [dir_validfids 'vnB\'];
dir_vnM = [dir_validfids 'vnM\'];

% list of included files (ref mk_validfids.m)
% NOTE: some files were dropped during manual cleaning! 
tmp = load([dir_validfids 'adjusted_filenames.mat']);
filenames = tmp.filenames;
clear tmp

% add path to some custom local functions
addpath([root 'code\1_BVP\localfunctions']);


%% Some settings from ex05_main

fs_org = 128;       % [Hz]
fs_int = nan;       % [Hz]

if isnan(fs_int)
    fs = fs_org;
else
    fs = fs_int;
end

fixed_time = 300;   % [s] Default 300s
fixed_size = fixed_time*fs;


%% Quick temporary inspection

% set true if wanted, skip otherwise
% .tif files are much more enjoyable to work with...
quick_inspection = false; 

if quick_inspection
    figure
    for i = 1:length(filenames)

        % get annotations
        tmpA = load([dir_vnA filenames{i} '_vnA.mat']);
        tmpB = load([dir_vnB filenames{i} '_vnB.mat']);
        tmpM = load([dir_vnM filenames{i} '_vnM.mat']);
        vnA = tmpA.vnA;
        vnB = tmpB.vnB;
        vnM = tmpM.vnM;

        % get PPG signal
        ppg = get_ppg(dir_original,filenames{i},fs_org,fs_int,fixed_size);

        hold on
        plot(ppg)
        plot(vnA,ppg(vnA),'r*')
        plot(vnB,ppg(vnB),'g*')
        plot(vnM,ppg(vnM),'b*')
        title([filenames{i} ' - latest annotations'])
        legend('nA','nB','nM')
        hold off

        % hold until done
        gotonext = false;
        while gotonext==false
            x = input('Continue? y/n [no]: ','s');
            if strcmp(x,'y')
                gotonext=true;
            end
        end

        % clear last figure
        clf;
    end
end


%% Save selection as matlab figures manually
% use breakpoint in for loop to cycle through selection

selection = filenames(1); % just use first for now
target_path = [dir_validfids 'mat_validfids\'];

figure
for i = 1:length(selection)
    
    printable_filename = strrep(selection{i},'_',' ');
    
    % get annotations
    tmpA = load([dir_vnA selection{i} '_vnA.mat']);
    tmpB = load([dir_vnB selection{i} '_vnB.mat']);
    tmpM = load([dir_vnM selection{i} '_vnM.mat']);
    vnA = tmpA.vnA;
    vnB = tmpB.vnB;
    vnM = tmpM.vnM;

    % get PPG signal
    ppg = get_ppg(dir_original,filenames{i},fs_org,fs_int,fixed_size);

    hold on
    plot(ppg)
    plot(vnA,ppg(vnA),'r*')
    plot(vnB,ppg(vnB),'g*')
    plot(vnM,ppg(vnM),'b*')
    title([printable_filename ' - full annotations'])
    xlabel('time [samples @128Hz]')
    ylabel('amplitude [a.u.]')
    legend('PPG','nA','nB','nM')
    hold off
end


%% Save plots as .tif for later inspection

% specify where to save .tif files
dir_tif = [dir_validfids 'tif_validfids\'];

for i = 1:length(filenames)
    
    printable_filename = strrep(filenames{i},'_',' ');
    disp(printable_filename)
    
    % get annotations
    tmpA = load([dir_vnA filenames{i} '_vnA.mat']);
    tmpB = load([dir_vnB filenames{i} '_vnB.mat']);
    tmpM = load([dir_vnM filenames{i} '_vnM.mat']);
    vnA = tmpA.vnA;
    vnB = tmpB.vnB;
    vnM = tmpM.vnM;

    % get PPG signal
    ppg = get_ppg(dir_original,filenames{i},fs_org,fs_int,fixed_size);
    
    % build figure
    figure('visible','off','name',['Cleaned BVP annotations for ' printable_filename]); 
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 50 25], 'PaperUnits', 'centimeters', 'PaperSize', [50 25]);

    hold on
    plot(ppg)
    plot(vnA,ppg(vnA),'r*')
    plot(vnB,ppg(vnB),'g*')
    plot(vnM,ppg(vnM),'b*')
    title([printable_filename ' - latest annotations'])
    legend('BVP','nA','nB','nM')
    hold off

%     % save figure as .tif
%     filepath = [dir_tif filenames{i}];
%     print(gcf,filepath,'-dtiff','-r300');
end


%% temp figs for rdeco_out (rnA) annotations

% specify where to save .tif files
dir_tif = [dir_validfids 'tif_rdeco_out\'];

for i = 1:length(filenames)
    
    printable_filename = strrep(filenames{i},'_',' ');
    disp(printable_filename)
    
    % get annotations
    tmp = load([dir_rdeco_out filenames{i} '_nA.mat']);
    rnA_sec = seconds(tmp.data.R_loc{:});   % [seconds]
    rnA = round(rnA_sec*fs);                % [samples]
    
    % get PPG signal
    ppg = get_ppg(dir_original,filenames{i},fs_org,fs_int,fixed_size);
    
    % build figure
    figure('visible','off','name',['Cleaned BVP annotations for ' printable_filename]); 
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 50 25], 'PaperUnits', 'centimeters', 'PaperSize', [50 25]);

    hold on
    plot(ppg)
    plot(rnA,ppg(rnA),'r*')
    title([printable_filename ' - manual annotations (import)'])
    legend('BVP','nA')
    hold off

%     % save figure as .tif
%     filepath = [dir_tif filenames{i}];
%     print(gcf,filepath,'-dtiff','-r300');
end


