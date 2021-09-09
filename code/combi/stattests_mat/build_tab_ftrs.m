close all
clear

%% Feature analysis - build feature list for appendix


%% dirs

root = 'C:\Users\jefma\repos\mox_thesis\';
datadir = [root 'data\'];

% save/load intermediate results in data folder
dir_intermediates = [datadir 'combi\intermediates_analysis\'];

% utility functions
% localfunctions self-defined within repo
% https://github.com/bastibe/Violinplot-Matlab
% https://github.com/eliduenisch/latexTable
addpath([root 'code\combi\stattests_mat\localfunctions\']);
addpath('C:\Users\jefma\repos\Violinplot-Matlab');
addpath('C:\Users\jefma\repos\latexTable');


%% Load data from previous steps and prepare for testing formats

% load relevant source data
orglog      = load([datadir 'logdata_all.mat']); % original ftr data
placebo     = load([datadir 'placebo.mat']);
attachment  = load([datadir 'attachment.mat']).attachment;
filenames   = load([datadir 'filenames.mat']).filenames;

% load prepped data from build_tfdata.m
SOURCE = load([dir_intermediates 'normtf_data.mat']);
tf_info = SOURCE.transformations;
ftr_ref = SOURCE.ftr_ref;
DATA = SOURCE.TFDATA;


%% 

printable = strrep(ftr_ref.names','_',' ');
feature = pad(printable);
ID = ftr_ref.ID_str';

for i=1:length(ID)
    if contains(feature{i},'PRV')
        disp([ID{i} ' & ' feature{i} ' &  \\'])
    end
end

fprintf('\n');

for i=1:length(ID)
    if contains(feature{i},'CDA')
        disp([ID{i} ' & ' feature{i} ' &  \\'])
    end
end

fprintf('\n');

for i=1:length(ID)
    if contains(feature{i},'DDA')
        disp([ID{i} ' & ' feature{i} ' &  \\'])
    end
end

fprintf('\n');

for i=1:length(ID)
    if contains(feature{i},'RSP')
        disp([ID{i} ' & ' feature{i} ' &  \\'])
    end
end