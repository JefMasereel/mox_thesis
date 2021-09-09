function handlers = build_handlers(log,placebo,attachment,filenames)

% log    struct  input data representing ftr values with headers
% dftr   struct  treatment effect data built from ftr data

%% Easier indexing formats
% pids   struct  subject ID indexing
% idx    struct  basic data indexing
% aidx   struct  SAAM scores indexing
% didx   struct  dftr indexing (treatment effects pre-post)
% ridx   struct  bcsr indexing (stress response rs-eye)

% local util functions to build dftr
root = 'C:\Users\jefma\repos\mox_thesis\';
addpath([root 'code\combi\stattests_mat\localfunctions\']);


%%  Standard feature values (ftr)

disp(placebo.labels)

pids.full     = placebo.labels.pid;
pids.oxytocin = placebo.labels.pid(contains(placebo.labels.treatment,'oxytocin'));
pids.placebo  = placebo.labels.pid(contains(placebo.labels.treatment,'placebo'));
idx.oxytocin  = contains(filenames,pids.oxytocin);
idx.placebo   = contains(filenames,pids.placebo);
idx.RS   = contains(filenames,'_RS_');
idx.EYE  = contains(filenames,'_EYE_');
idx.PRE  = contains(filenames,'PRE');
idx.POST = contains(filenames,'POST');


%% Treatment effect data (dftr)
% dftr = ftr(post) - ftr(pre)

% dftr      treatment effect data (simple form)
% dftr_ext  extended form, includes condition, pid & placebo columns
dftr     = build_dftr(log,filenames,placebo);
dftr_ext = build_dftr_ext(log,attachment,filenames,placebo);

% indexing formats for treatment effect data
% note: transposed to obtain column arrays (matches ftr formatting)
didx.oxytocin = contains(dftr.rownames,pids.oxytocin)';
didx.placebo = contains(dftr.rownames,pids.placebo)';
didx.RS = contains(dftr.rownames,'_RS')';
didx.EYE = contains(dftr.rownames,'_EYE')';


%% Baseline-corrected stress response data (bcsr)
% bcsr = ftr(eye) - ftr(rest)

% bcsr = ftr(eye) - ftr(rest)
bcsr = build_bcsr(log,placebo);

% indexing formats
ridx.oxytocin = contains(bcsr.rownames,pids.oxytocin)';
ridx.placebo  = contains(bcsr.rownames,pids.placebo)';
ridx.PRE  = contains(bcsr.rownames,'_PRE')';
ridx.POST = contains(bcsr.rownames,'_POST')';


%% Treatment effect on bcsr (dbcsr)
% dbcsr = bcsr(post) - bcsr(pre)

% bcsr = ftr(eye) - ftr(rest)
dbcsr = build_dbcsr(bcsr,log,placebo);

% indexing formats
dridx.oxytocin = contains(dbcsr.rownames,pids.oxytocin)';
dridx.placebo  = contains(dbcsr.rownames,pids.placebo)';


%% Easier indexing formats for attachment data (missing pids!)

aidx.RS  = true(size(attachment.scores,1),1);
aidx.EYE = true(size(attachment.scores,1),1);
aidx.oxytocin = (attachment.scores(:,2)==1);
aidx.placebo = (attachment.scores(:,2)==2);

% Warning: missing pids for certain samples!
aidx.RS(40) = false;
aidx.EYE([1,7,37,38,40]) = false;

% % indexing attachment scores for treatment effect correlations
% daidx.RS  = true(1,size(attachment.scores,1));
% daidx.EYE = true(1,size(attachment.scores,1));
% daidx.RS(40) = [];
% daidx.EYE([1,7,37,38,40]) = [];


%% collect results, return as struct

handlers = struct('pids',pids,...
                  'idx',idx,'aidx',aidx,...
                  'didx',didx,'dftr',dftr,'dftr_ext',dftr_ext,...
                  'ridx',ridx,'bcsr',bcsr,...
                  'dridx',dridx,'dbcsr',dbcsr);

% when implementing in script, select relevant data:

% hdl = build_handlers(log,placebo,attachment,filenames);
% pids = hdl.pids;
% idx = hdl.idx;
% aidx = hdl.aidx;
% didx = hdl.didx;
% dftr = hdl.dftr;
% dftr_ext = hdl.dftr_ext;
% ridx = hdl.ridx;
% bcsr = hdl.bcsr;
% dridx = hdl.dridx;
% dbcsr = hdl.dbcsr;

end