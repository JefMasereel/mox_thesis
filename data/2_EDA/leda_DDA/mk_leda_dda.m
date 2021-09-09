close all
clear

%% EDA - Run DDA Batch analysis with Ledalab package

% using Ledalab package for EDA decomposition
% git clone https://github.com/ledalab/ledalab.git
% documentation and research papers on http://ledalab.de
addpath('C:\Users\jefma\repos\ledalab\');

% NOTE
% copy the contents of leda_in/ to leda_dda/
% then run a copy of this script inside leda_dda/
% easiest way to apply batch analysis with Ledalab
%   (Ledalab spits out result files in pwd...)


%% Apply Batch Analysis - DDA method
% Discrete Decomposition Analysis
% "This method is especially suited for studies of physiological models 
% of the SCR but may be slow for large data."

% % Args
% open            'mat' for new files, once processed they are default 'leda'
% analyze         use CDA as recommended, more interesting than DDA (and faster)
% optimize        3 iterations should be enough to get an acceptable error
% export_era      screen opens 5s at onset event, recommended 1-4s after impulse event
%                 --> choose response window 1-10 seconds?
%                 amp_threshold default .01, looks alright
% export_scrlist  use same settings as export_era, should be fine
% overview        set 1 to save overview images of each file

% batch analysis for all files in working directory
Ledalab(pwd,'open','mat','analyze','DDA','optimize',3,...
        'export_era',[1 4 .01 1], 'export_scrlist',[.05 1],...
        'overview',1);


%% Apply EDA_Results.m
% note: in automation easier to recalculate scores in script
% ref EDA_Results.m, simple for loop can be reused and extended

% EDA_Results; % from cmd line, in target pwd


