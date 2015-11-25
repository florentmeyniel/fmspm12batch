% script to 'custumize' the spm_check_results.m function from SPM
%   -> request a rfx model folder (to get the model name)
%   -> request a list of suject
%   -> request which estimated contrast to plot in the model
%   -> request how to threshold the map
% It then finds the individual SPM.mat and plot the Maximum Intensity
% Projections, for all subjects.

% initialize the SPM window.
try ver = spm('Version'); 
catch
    error('the SPM toolbox is not in the Matlab''s path')
end
spm_results_ui('clear')
spm_results_ui('close')

% get the rfx model directory (to get the model name)
modeldir = spm_select(1,'dir','Select rfx model');
ind = strfind(modeldir, '/');
modelname = modeldir(ind(end)+1:end);
ind = strfind(modeldir, '/group_analysis/');
rootdir = modeldir(1:ind-1);

% get list of subject to plot
sublist = spm_input('subject list','+1','e','');

% built list of individual SPM.mat
SPMs = [];
for iSub = 1:numel(sublist)
    SPMs(iSub,:) = ...
        sprintf('%s/subj%02.0f/first_level_estimates/%s/SPM.mat', ...
        rootdir, sublist(iSub), modelname);
end
SPMs = char(SPMs);

% check that files exist
for iSub = 1:numel(sublist)
    if ~exist(SPMs(iSub,:), 'file')
        error('the following file does not exist: \n %s', SPMs(iSub,:))
    end
end

% Make the multiple plot.
spm_check_results(SPMs, [])