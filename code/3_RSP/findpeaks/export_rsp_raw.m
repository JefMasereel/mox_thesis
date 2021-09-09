close all
clear


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\data\';
dir_rsp = [root '3_RSP\findpeaks\'];
dir_from = [root 'ORIGINAL\rawdata_signals\'];
dir_to = [dir_rsp 'export_rsp\raw\'];

% list of included files (ref init_process.m)
% NOTE: some files were dropped during manual cleaning!
tmp = load([dir_rsp 'adjusted_filenames.mat']);
filenames = tmp.filenames;
clear tmp


%% Iterate over source files and save to .mat

fs = 128;               % [Hz]
fixlength = 300*fs;     % [s]

for i = 1:length(filenames)
    
    printable_filename = strrep(filenames{i},'_',' ');
    disp(printable_filename)
    
    tmp = readtable([dir_from filenames{i} '.txt']);
    rsp_raw = tmp.Var4';
    
    rsp = minimal_prep(rsp_raw,fixlength);

    savepath = [dir_to filenames{i} '.mat'];
    save(savepath,'rsp','fs');
end


%% Local function

function rsp = minimal_prep(rsp,fixlength)
    
    % Set fixed length @ 300 seconds
    try
        rsp = rsp(1:fixlength);
    catch
        fprintf('signal has only %d samples \n',length(rsp));
    end

    % remove nan values at end
    rsp = rsp(~isnan(rsp));  
end