% EXAMPLE
% Script to specify the design matrix (experiment-specific and subject-specific) for the 
% first level analysis of fMRI data.
%
% This script allows the specification of parametric modulations. If there is no parametric 
% modulation of a given onset regressor (e.g. the 2nd one), just leave the fields empty:
%         pmod(2).name{1}     = [];
%         pmod(2).param{1}    = [];
%         pmod(2).poly{1}     = [];
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
addpath ../
tmp = mfilename;
modelname = tmp(1:end-10);

% Parameters to specify. NB: session 1 is the training session
fmrisess = 2:5;

for iSub = 1:length(sublist)
    fprintf('\n Computing subject %d (ID: subj %02.0f)', iSub, sublist(iSub))
    for iSess = 1:length(fmrisess);
        
        clear names onsets durations pmod
        
        % initialize pmod structure
        pmod = struct('name',{''},'param',{},'poly',{});
        
	% load data for this subject.
        opt.datadir = behav_dir;
        opt.sublist = {sprintf('s%02.0f', sublist(iSub))};
        opt.sesslist = fmrisess(iSess);
        out = m20151119_getData_v2(opt);
        
        
        % Reg1: Stimulus Onset
        % =================================================================
        names{1}            = 'Stimulus Onset';
        
        ons_StimOn = out.t_AnimOnset{1} - out.t_T0{1};
        
        onsets{1}            = ons_StimOn;
        durations{1}         = out.t_AnimOffset{1} - out.t_AnimOnset{1};
        
        Left_lum  = out.Left_lum{1};
        Right_lum = out.Right_lum{1};
        
        pmod(1).name{1}     = 'Right Luminance';
        pmod(1).param{1}    = zscore(Right_lum);
        pmod(1).poly{1}     = 1;
        pmod(1).name{2}     = 'Left Luminance';
        pmod(1).param{2}    = zscore(Left_lum);
        pmod(1).poly{2}     = 1;
        
        
        % Reg2: Response Onset
        % =================================================================
        names{2}            = 'Response Onset';
        
        ons_RepConf = out.t_Response{1} - out.t_T0{1};
        isValid = ~isnan(ons_RepConf);
        
        onsets{2}           = ons_RepConf(isValid);
        durations{2}        = 0;
        pmod(2).name{1}     = 'RT';
        pmod(2).param{1}    = zscore(out.RT{1}(isValid));
        pmod(2).poly{1}     = 1;
        pmod(2).name{2}     = 'Right > Left answer';
        pmod(2).param{2}    = out.Choice{1}(isValid);
        pmod(2).poly{2}     = 1;
        
        % Reg3: Feedback Onset
        % =================================================================
        names{3}            = 'Feedback Onset';
        
        ons_FB = out.t_feedback{1} - out.t_T0{1};
                
        onsets{3}           = ons_FB(isValid);
        durations{3}        = 0;
        pmod(3).name{1}     = 'valence';
        pmod(3).param{1}    = out.Correct{1}(isValid);
        pmod(3).poly{1}     = 1;
                
                
        % Save multicond
        % =================================================================
        % Check that the pmod structure array is 1xn is size
        if length(pmod) ~= length(names)
            error('pmod is size %d whereas names is size %d', ...
                length(pmod), length(names))
        end
        
        fname = sprintf('%s/subj%02.0f/MultiCond/%s_multicond_session%d.mat', ...
            irm_dir, sublist(iSub), modelname, iSess);
        
        if ~exist(sprintf('%s/subj%02.0f/MultiCond/', irm_dir, sublist(iSub)))
            mkdir(sprintf('%s/subj%02.0f/MultiCond/', irm_dir, sublist(iSub)))
        end
        save(fname, 'names', 'onsets', 'durations', 'pmod')
    end
end


