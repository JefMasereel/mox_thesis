close all
clear

%% EDA - Prepare raw signal for processing with Ledalab

% using Ledalab package for EDA decomposition
% git clone https://github.com/ledalab/ledalab.git
% documentation and research papers on http://ledalab.de


%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\data\';
dir_org = [root 'ORIGINAL\rawdata_signals\'];
dir_eda = [root '2_EDA\'];

% import toolkit functions
addpath('C:\Users\jefma\repos\ledalab\');    % Ledalab repo

% files to be included
% NOTE: run init_process.m to build headers and filenames
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
save([dir_eda 'adjusted_filenames.mat'],'filenames');


%% Read data and reformat to leda_in

% downsample
fs_org = 128; % [Hz]
fs_eda = 32;

% fix signal length
fixed_duration_RS = 300;
fixed_duration_EYE = 300;

% define labels of interest
closed_label = "RS232 Trigger: 67('C')";
closed_nid   = 1;
closed_name  = 'closed';
open_label   = "RS232 Trigger: 79('O')";
open_nid     = 2;
open_name    = 'open';

for i = 1:length(filenames)
    
    disp(filenames{i})
    
    tmp = readtable([dir_org filenames{i} '.txt']);
    EDA = tmp.Var2;
    
    % get event data
    EVTS = double.empty;
    if size(tmp,2)==9
        EVTS = tmp.Var9;
    end
    
    % check & clean nan-values
    if any(isnan(EDA))
        if ~any(find(isnan(EDA))<length(EDA)-1)
            EDA = EDA(~isnan(EDA));
            try EVTS = EVTS(~isnan(EDA)); end
        else
            disp('Measurement contains nan-values in body')
        end
    end
    
    % set standard signal lengths as wanted
    if contains(filenames{i},'_RS_')
        [EDA,~] = fixlength_eda(EDA,fs_org,fixed_duration_RS);
    elseif contains(filenames{i},'_EYE_')
        [EDA,~] = fixlength_eda(EDA,fs_org,fixed_duration_EYE);
    end
    
    % downsample all data to fs_eda (assuming fs_eda < fs_org)
    eda_idx     = 1:fs_org/fs_eda:length(EDA); % time@fs_org!
    eda_signal  = EDA(eda_idx);
    eda_time    = 1/fs_eda:1/fs_eda:length(eda_signal)/fs_eda;
    
    % rebuild data in new format
    data.conductance    = eda_signal';
    data.time           = eda_time;
    data.timeoff        = 0;
    data.event          = [];
    
    % rebuild event data
    for j = 1:length(EVTS)
        
        if isempty(EVTS{j})
            event_nid = nan;
        elseif EVTS{j}==closed_label
            event_nid  = closed_nid;
            event_name = closed_name;
        elseif EVTS{j}==open_label
            event_nid  = open_nid;
            event_name = open_name;
        else
%             disp('invalid event label found')
%             disp(EVTS{j})
            event_nid = nan;
        end
        
        if ~isnan(event_nid)
            
            % get time [sec] of event
            [~,kk] = min(abs(eda_idx-j));
            event_time = eda_idx(kk)/fs_org;
            
            % add event to data
            event_struct = struct('time',event_time,...
                                  'nid',event_nid,...
                                  'name',event_name,...
                                  'userdata',char.empty);
            data.event = [data.event event_struct];
        end        
    end
    
    % Save results
    save([dir_eda 'leda_in\' filenames{i} '.mat'],'data');
    
end


%% Local function

function [eda,flag_short] = fixlength_eda(eda,fs,fixed_duration)
% eda               time signal (electrodermal activity)
% fs                sample frequency of eda signal
% fixed_duration    standardize duration [sec] of measurements
% eda               EDA signal of fixed duration @fs
% flag_short        raise flag (bool) if EDA is shorter than specified size

flag_short = false;
try
    eda = eda(1:fixed_duration*fs);
catch
    flag_short = true;
    % return eda as is
end
end


