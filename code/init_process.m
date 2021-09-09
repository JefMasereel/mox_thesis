close all
clear

%% General initialization of signal processing

% Script based on previous efforts in qssr_code repo
% as part of thesis by Jef Masereel

% initial build of all overarching files, nl.
%   headers             overview of available files in dataset
%   filenames           easier iterable format to use in for loops
%   placebo_rawdata     
%   placebo_labels

% note: placebo_rawdata.txt has to be me


%% dir

% define project root and add to path for later use
root = 'C:\Users\jefma\repos\mox_thesis\';
addpath(root);

% directories for this script
dir_src = [root 'data\ORIGINAL\'];
dir_hdr = [root 'data\'];


%% specify and build iterable files

% list of measurements to be excluded due to artefacts [UNUSED]
% ignore_full = {'PP05_RS_PRE';'PP05_RS_POST';'PP05_EYE_PRE';'PP05_EYE_POST';...
%                'PP14_RS_PRE';'PP14_RS_POST';'PP14_EYE_PRE';'PP14_EYE_POST';...
%                'PP53_RS_PRE';'PP53_RS_POST';'PP53_EYE_PRE';'PP53_EYE_POST';...
%                'PP49_RS_PRE';'PP49_RS_POST';'PP49_EYE_PRE';'PP49_EYE_POST1';'PP49_EYE_POST2';...
%                'PP52_RS_PRE';'PP52_RS_POST';'PP52_EYE_PRE';'PP52_EYE_POST'};
% ignore doubles
ignore_doubles = {'PP31_EYE_POST2';
                  'PP32_EYE_PRE2';'PP32_EYE_POST2';
                  'PP33_EYE_PRE2';'PP33_EYE_POST2';
                  'PP34_EYE_PRE2';'PP34_EYE_POST2';
                  'PP35_EYE_PRE2';'PP35_EYE_POST2';
                  'PP36_EYE_PRE2';'PP36_EYE_POST2';
                  'PP38_EYE_POST2';
                  'PP48_EYE_PRE2';
                  'PP49_EYE_POST2';
                  'PP51_EYE_PRE2';
                  'PP54_EYE_POST2';
                  'PP55_EYE_POST2';
                  'PP56_EYE_PRE2'};
% ignore files with state=EYE
ignore_eye_only = {'PP01_EYE_PRE';'PP01_EYE_POST';
                   'PP07_EYE_PRE'; %'PP07_EYE_POST';
%                    'PP38_EYE_PRE';'PP38_EYE_POST'}
                   'PP37_EYE_PRE';'PP37_EYE_POST'};

% All relevant ignorelists concatenated
ignore_msrs = [ignore_doubles; ignore_eye_only];

% Start building iterables from there
headers     = init_headers([dir_hdr 'headers.mat']);
filenames   = init_filenames(headers,[dir_hdr 'filenames.mat'],ignore_msrs);
adj_files   = adjust_filenames(filenames,dir_src,dir_hdr);
placebo     = init_placebolabels(dir_src,dir_hdr);
attachment  = init_attachmentlabels(dir_src,dir_hdr);


%% Local functions

function headers = init_headers(savefile)
% Construct a new iterable structure to read out data
% lots of manual tweaking for irregularities in this dataset
% complex structure, prefer filenames.mat for simple iterations
% savefile = PATH to save results to, incl. filename and type
% f.e. '~/data/headers_dmy.mat'

    PREFIX = 'PP';

    % missing measurement indices        %% only 1 full pid missing
    missing_idx = [40];

    % available measurements
    msr_list = 1:56;
    msr_list(missing_idx) = [];
    tot_nb = length(msr_list);

    % num2string patient ID
    MSR_IDS = string(zeros(size(msr_list)));
    for i = 1:length(msr_list)
        MSR_IDS(i) = num2str(msr_list(i),'%02d');
    end

    % define label idx
    STATE = {'_RS_', '_EYE_'};
    TIME = {'PRE', 'POST'};

    % EYE data contains double measurements
    EYE_NB = ones(2,56);
    EYE_NB(1,[32:36 48 51 56]) = 2; % PRE
    EYE_NB(2,[31:36 38 49 54 55]) = 2; % POST
    EYE_NB(:,missing_idx) = [];
    SUFFIX = {'1', '2'};

    HEADERS = cell(2,tot_nb,2); % (state,msr_id,time,doubles)
    % iterate over all elements
    for state = 1:2                          % rest, eye
        for id = 1:length(MSR_IDS)           % 1:56 minus missing data
            for time = 1:2                   % pre, post
                for nb = 1:2                 % double measurements present

                    if EYE_NB(time,id)==2 && state==2
                        header = join([PREFIX MSR_IDS{id} STATE{state} TIME{time} SUFFIX{nb}],'');
                    elseif (state==1 || EYE_NB(time,id)==1) && nb==2
                        header = nan;
                    else
                        header = join([PREFIX MSR_IDS{id} STATE{state} TIME{time}],'');
                    end

                    HEADERS{state,id,time,nb} = header;

                end
            end
        end
    end
    
    % results
    headers = HEADERS;
    save(savefile,'HEADERS')
end


function filenames = init_filenames(headers,savefile,ignore_msrs)
% requires headers file
% make that with init_headers.m

% hdr_file      path of header file, should include .mat suffix
% savefile      path to save results to, incl. filename and type
% ignore_msrs   list of filenames to be removed from further processing

    nb_stts  = size(headers,1); % STATE (rest, eye)
    nb_ids   = size(headers,2); % ID
    nb_times = size(headers,3); % TIME (pre post)
    nb_dbls  = size(headers,4); % double measurements
    
    filenames = {};
    for state = 1:nb_stts
        for id = 1:nb_ids
            for time = 1:nb_times
                for nb = 1:nb_dbls

                    header = headers{state,id,time,nb};

                    if not(isnan(header))                    
                        filenames{end+1,1} = header;

                    end
                end
            end
        end
    end

    % exclude specified measurements (if ignore_msrs not empty)
    filenames = setdiff(filenames, ignore_msrs);

    % save results
    save(savefile,'filenames')
end


function filenames = adjust_filenames(filenames, dir_src, dir_hdr)
% check for missing filenames and correct iterable

    adj_files = {};
    for i = 1:length(filenames)
        try 
            tmp = readtable([dir_src 'rawdata_signals\' filenames{i} '.txt']);
            adj_files{end+1} = filenames{i};
        catch
            % if file is missing, filename is not added to adj_files list
        end
    end
    fprintf("%d missing files in dataset, removed from further analysis \n",...
        length(filenames)-length(adj_files));
    disp(setdiff(filenames,adj_files))
    
    % save results
    filenames = adj_files'; % transposed for formatting
    save([dir_hdr 'adjusted_filenames.mat'],'filenames');
end

function placebo = init_placebolabels(dir_src,dir_hdr)
% dir_rawfile   bin placebo labels (.txt copied from excel overview)
% RESULT        saves labels struct with pids and matching placebo label

    % load source information (might require manual prep)
    raw = load([dir_src 'rawdata_placebo.txt']);

    % reformatted data structure
    nb_patients = size(raw,1);
    labels.pid = cell(nb_patients,1);       % patient identifiers
    labels.bin = zeros(nb_patients,1);      % binary placebo labels
    labels.treatment = cell(nb_patients,1); % label interpretation

    % link to patient ID
    for id = 1:nb_patients
        nb = num2str(id,'%02d');
        labels.pid{id} = ['PP' nb];
        labels.bin(id) = raw(id);
        if raw(id)==1
            labels.treatment{id} = 'oxytocin';
        else
            labels.treatment{id} = 'placebo';
        end
    end
    
    % FOR DEMO ONLY (reduced dataset to single participant)
    labels.pid = cell(56,1);
    labels.bin = cell(56,1); 
    labels.treatment = cell(56,1);
    labels.pid{2} = 'PP02';
    labels.bin{2} = 0;
    labels.treatment{2} = 'oxytocin';
    
    % collect results
    placebo = labels;
    save([dir_hdr 'placebo.mat'],'labels');
end


function attachment = init_attachmentlabels(dir_src,dir_hdr)

    % load source information (might require manual prep)
    SAAM_raw = importdata([dir_src 'rawdata_attachment.txt']);
    
    % reformat SAAM score data
    % formatting similar to placebo for consistency
    attachment.names = {'pid','placebo','anxiety','security','avoidance'};
    attachment.scores = zeros(56,5);

    for i = 1:length(SAAM_raw)/5

        pid         = SAAM_raw(5*i-4); % subject ID
        spray       = SAAM_raw(5*i-3); % oxytocin = 1, placebo = 2
        anxiety     = SAAM_raw(5*i-2); % ordinal scale 1-7
        security    = SAAM_raw(5*i-1); % ordinal scale 1-7
        avoidance   = SAAM_raw(5*i);   % ordinal scale 1-7

%         assert(i==pid); % commented for demo with devalued data!

        attachment.scores(pid,1) = int8(pid);
        attachment.scores(pid,2) = int8(spray);
        attachment.scores(pid,3) = anxiety;
        attachment.scores(pid,4) = security;
        attachment.scores(pid,5) = avoidance;
    end
    
    savepath = [dir_hdr 'attachment.mat'];
    save(savepath,'attachment');
end