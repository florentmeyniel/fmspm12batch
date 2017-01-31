% EXAMPLE
% Script to specify the contrast (experiment-specific, subject-specific and model-specific) 
% for the first level analysis of fMRI data.
%
% Note that the file prefix should be matched: 
% 	prefix_contrast.m 
%	prefix_multicond.m 
% 	prefix_paramfile.m

% Specify directory to use
behav_dir   = '/neurospin/unicog/protocols/IRMf/Meyniel_LumiConfidence_2014/Behavioral_data/';
irm_dir     = '/neurospin/unicog/protocols/IRMf/Meyniel_LumiConfidence_2014/MRI_data/analyzed_data/';
spm_dir     = '/volatile/meyniel/toolbox/matlab/spm12/';

% Initialize script
addpath(spm_dir)
locmodelname = mfilename;
locmodelname = locmodelname(1:end-9);

% Parameters to specify
fmrisess = 2:5;
nsess    = numel(fmrisess);

% initialize
clear names type values

% initialize contrast
names{1}  = 'Stim ons.';
names{2}  = 'Stim ons. x Right Lum';
names{3}  = 'Stim ons. x Left Lum';
names{4}  = 'Stim ons. x R > L Lum';
names{5}  = 'Stim ons. x R < L Lum';
names{6}  = 'Motor response ons.';
names{7}  = 'Motor response ons. x RT';
names{8}  = 'Motor response ons. x (-1) RT';
names{9}  = 'Motor response ons. x R > L';
names{10} = 'Motor response ons. x R < L';
names{11} = 'FB ons. x pos > neg';
names{12} = 'FB ons. x pos < neg';

for iCon = 1:numel(names)
    type{iCon} = 'T';
end

% GET PADDING FOR CONFOUND REGRESSORS
% for debug mode: call model parameters. In the standard use, these
% parameters are already loaded in the workspace when this script is
% called.
% eval(sprintf('%s_paramfile', locmodelname));

if physiocorr.include == 0
    % include the movement regressors
    padding = zeros(1, 6);
else
    % add sin and cos wave modelling physiological regressors
    n_TAPAS = ...
        2*physiocorr.opt.order_cardiac + ...
        2*physiocorr.opt.order_resp + ...
        4*physiocorr.opt.order_interaction;
    
    % add the other regressors (if any)
    switch physiocorr.type
        case 'RETROICOR_HR'
            n_TAPAS = n_TAPAS + 1;
        case 'RETROICOR_RV'
            n_TAPAS = n_TAPAS + 1;
        case 'RETROICOR_HR_RVT'
            n_TAPAS = n_TAPAS + 2;
    end
    
    % include the movement + physiological regressors
    padding = zeros(1, 6+n_TAPAS);
end

for iSub = 1:length(sublist)
    
    fprintf('\n Computing subject %d', iSub)
    clear values
    for iCon = 1:numel(names)
        values{iCon} = [];
    end
        
    
    for iSess = 1:nsess
        
        % Stim ons.
        values{1} = [
            values{1}, ...
            1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0, ...
            padding];               % confounds
        
        % 'Stim ons. x Right Lum';
        values{2} = [
            values{2}, ...
            0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0, ...
            padding];               % confounds
        
        % 'Stim ons. x Left Lum';
        values{3} = [
            values{3}, ...
            0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0, ...
            padding];               % confounds
        
        % 'Stim ons. x R > L Lum';
        values{4} = [
            values{4}, ...
            0 0 1 0 -1 0 0 0 0 0 0 0 0 0 0 0, ...
            padding];               % confounds
        
        % 'Motor response ons.';
        values{6} = [
            values{6}, ...
            0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0, ...
            padding];               % confounds
        
        % 'Motor response ons. x RT'
        values{7} = [
            values{7}, ...
            0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0, ...
            padding];               % confounds
        
        % 'Motor response ons. x R > L';
        values{9} = [
            values{9}, ...
            0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0, ...
            padding];               % confounds
        
        % 'FB ons. x pos > neg';
        values{11} = [
            values{11}, ...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0, ...
            padding];               % confounds
        
    end
    
    % Compute the reverse contrasts
    values{5}  = -values{4};
    values{8}  = -values{7};
    values{10}  = -values{9};
    values{12} = -values{11};
    
    % Normalize contrasts
    for iCon = 1:numel(names)
        tmp = values{iCon}(values{iCon} ~= 0);
        
        if all(tmp>0)
            % ensure that sum(c) = 1
            values{iCon} = values{iCon} / sum(tmp);
        elseif all(tmp<0)
            % ensure that sum(c) = -1
            values{iCon} = values{iCon} / sum(abs(tmp));
        else
            % ensure that sum(abs(c)) = 1 and sum(c) = 0
            npos = sum(tmp>0);
            nneg = sum(tmp<0);
            
            wneg = 1 / (2*nneg);
            wpos = wneg * nneg/npos;
            
            tmp(tmp>0) = tmp(tmp>0) * wpos;
            tmp(tmp<0) = tmp(tmp<0) * wneg;
            
            values{iCon}(values{iCon} ~= 0) = tmp;
        end
    end
    
    % append zeros for the session constant
    for iCon = 1:numel(names)
        values{iCon} = [values{iCon}, zeros(1, numel(fmrisess))];
    end
    
    % Save contrast
    % =====================================================================
    fname = sprintf('%s/subj%02.0f/MultiCond/%s_contrast.mat', ...
        irm_dir, sublist(iSub), locmodelname);
    
    if ~exist(sprintf('%s/subj%02.0f/MultiCond/', irm_dir, sublist(iSub)))
        mkdir(sprintf('%s/subj%02.0f/MultiCond/', irm_dir, sublist(iSub)))
    end
    save(fname, 'names', 'type', 'values')
end
fprintf('\n')

