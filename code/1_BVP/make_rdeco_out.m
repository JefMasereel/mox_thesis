close all
clear

%% BVP - tool to speed up manual cleaning in RDECO

% ref RDECO by Jonathan Moeyersons:
% https://physionet.org/content/r-deco/1.0.0/

% NOTE
% works best when running a copy of this script inside the rdeco_out folder
% script loads files from rdeco_in --> save some clicks every cycle


%% dir

root = 'C:\Users\jefma\repos\mox_thesis\';

dir_hdr = [root 'data\1_BVP\'];
dir_rdeco_in  = [dir_hdr 'rdeco_in\'];
dir_rdeco_out = [dir_hdr 'rdeco_out\'];

% git clone https://gitlab.esat.kuleuven.be/biomed-public/r-deco.git
addpath('C:\Users\jefma\repos\r-deco')

% list of included files
tmp = load([dir_hdr 'adjusted_filenames.mat']);
filenames = tmp.filenames;
clear tmp


%% Use a progress file for intermittent sessions (optional)

% uncomment next line to refresh the savefile
init_progresstable(filenames,pwd)

% choose fiducial to work on
fiducial = 'nA';

% load progress from last saved session
progress = readtable([pwd '\cleaning_progress.txt']);
start_i = find(progress{:,fiducial}==0,1,'first');


%% Start cleaning

R_DECO

if isempty(start_i)
    disp('All files have been cleaned for this fiducial')
else
    for i = start_i:size(progress,1)
        
        % load file and apply R_DECO
        header = progress{i,'filenames'}{1};
        disp(header)
        clipboard('copy',[header '_' fiducial]);

        tmp = load([dir_rdeco_in header '.mat'],'ppg');
        ppg = tmp.ppg;
        clear tmp
        
        % update progress file (temp, needs to be reviewed and saved)
        progress{i,fiducial} = 1;
        
    end
end

disp('full column done, confirm and save table')
disp(progress)

prompt = 'Content ok? Y/N [N]: ';
str = input(prompt,'s');
if str=='Y'
    writetable(progress,[pwd '\cleaning_progress.txt'],'Delimiter','\t');
end


%% manual save

% progress{1:100,'nA'} = ones(100,1);
% writetable(progress,[pwd 'cleaning_progress.txt'],'Delimiter','\t');


%% Local functions

function init_progresstable(filenames,dir)
% set up a fresh .txt table to keep track of manual cleaning progress
% can be helpful when cleaning over interrupted sessions

% specify desired columns
i  = transpose(1:length(filenames));
nA = zeros(length(filenames),1);
mid = filenames';

progress = table(i,mid, nA);
writetable(progress,[dir '\cleaning_progress.txt'],'Delimiter','\t');
end


