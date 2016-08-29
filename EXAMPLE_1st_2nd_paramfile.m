% Script to set the experiment-specific pamareters for the first and second
% level analysis of fMRI data.

% For the batch system to work properly, a specific folder organization,
% with particular names, is required.
% Starting from the root directory:
%   datadir/
%     |
%     |_ subj01/
%     |_ subj02/
%     |_ ...
%     |_ group_analysis/rfxmodel/
%
%   datadir/
%     |
%     |_ subj01/
%       |_ first_level_estimates/   => where batches, beta, con, spmT, ... are saved
%       |_ MutliCond/               => where regressors & contrast specification is saved
%       |_ preprocAnat/             => preprocessed anatomy (for visualization)
%       |_ preprocEPI/              => preprocessed EPIs, movement parameters and SliceTimingInfo.mat
%
%

% Specification of parameters
% =========================================================================

% Option for distributed computing
useparallel.do    = 1;               % 1 to run in parallel, 0 to run serially
useparallel.max   = 6;               % maximum number of batch launched in parallel
useparallel.cmd   = 'matlab'; % command to invoke matlab from a terminal

% subject list
sublist           = [1:11];

% Actions to perform
% Possible actions are: 'specify', 'estimate', 'contrast', 'smooth', 'rfx_ttest', 'run'
% Note that 'run' runs the specified batch (without 'run', the batch is simply saved)
% 'specify', 'estimate' and 'contrast' are implemented with fmspm12batch_1st_level
% 'smooth' and 'rfx_ttest' are implemented with fmspm12batch_2nd_level
actions           = {'rfx_ttest', 'run'};

% locate the data
spm_path          = '/volatile/meyniel/toolbox/matlab/spm12';
datadir           = '/neurospin/unicog/protocols/IRMf/Meyniel_MarkovGuess_2014/MRI_data/analyzed_data';
regexp_func       = '^sw.*\.nii';                     % regular expression to recognize functional sessions to analyze
funcdir           = 'preprocEPI';                     % name of the folder containing EPI


% ----- 1st level parameters -----

% temporal basis function to use.
% 'hrf', 'hrf+deriv', 'hrf+2derivs', 'fir', 'fourier'
bases_functions   = 'hrf+deriv';

% Delete previous contrast (1 for yes, 0 for no)
deleteprevcon     = 1;

% reference slice for microtime resolution. If a slice timing correction
% was applied, with 1st slice at 0ms, the refslice should be aligned on the
% first time bin (=1) of the microtime resolution.
refslice          = 1;

% Mask (leave empty quotes '' to have no masking)
% spm_path/tpm/mask_ICV.nii,1 is a normalized, in-brain mask.
mask              = sprintf('%s/tpm/mask_ICV.nii,1', spm_path);

% Inclusion of physiological artifacts as co-variates
% The cardiac and respiratory phases (and interactions) can be included at
% various orders: this is the standard RETROICOR correction. In addition,
% the heart rate variability (HR) and respiratory volume per time (RVT) may
% also be included. Regressors were used using TAPAS during the
% preprocessing step.
physiocorr.include               = 1; % 1 or 0 for yes or no
physiocorr.type                  = 'RETROICOR'; % 'RETROICOR', 'RETROICOR_HR', 'RETROICOR_RVT', 'RETROICOR_HR_RVT'
physiocorr.opt.order_cardiac     = 3; % order of the Fourier expansion for the cardiac phase
physiocorr.opt.order_resp        = 4; % order of the Fourier expansion for the respiratory phase
physiocorr.opt.order_interaction = 1;

% ----- 2nd level paramters -----

% List of contrast computed. The number of the contrasts are as they can be
% read in 1st level SPM.mat.
% To compute all contrasts, leave an empty vector []
con_list          = [1 3 4];

% smoothing kernel for the 2nd level
smoothing_kernel2 = [6 6 6];

% type of 1st level contrast should be taken into the 2nd level analysis
% 'con' or 'scon' for unsmoothed & smoothed contrast
contrast_type     = 'scon' ;

% Mask (leave empty quotes '' to have no masking) for 2nd level
% spm_path/tmp/mask_ICV.nii,1 is a normalized, in-brain mask.
% mask_GreyMatter.nii,1 is a normalized, grey matter mask (made with ImCalc
%   with spm_path/tmp/mask_TPM.nii,1 > 0.05
mask2             = sprintf('%s/group_analysis/template/mask_GreyMatter.nii,1', datadir);


% Parse parameters
% =========================================================================

% lower the case of actions to perform
actions = lower(actions);

% get model name automatically
fname     = mfilename;
modelname = fname(1:end-10); % assume that file name is in the form modelname_paramfile.m

% get the fMRI acquisition TR
nSub = length(sublist);
TR   = zeros(1, nSub);
for iSub = 1:nSub
    fname = sprintf('%s/subj%02.0f/%s/SliceTimingInfo.mat', datadir, sublist(iSub), funcdir);
    if ~exist(fname, 'file')
        error('could not find the mat file with acquisition parameters at %s', fname)
    else
        dat = load(fname);
    end
    TR(iSub) = dat.TR;
end

% Check wether all participants have the same TR
if length(unique(TR)) == 1
    % all participants have the same TR => set it in a 1-value vector
    TR = unique(TR);
end

% Check whether masks exist.
if ~exist(mask(1:end-2))
    error('the mask file does not exist: %s', mask(1:end-2))
end
if ~exist(mask2(1:end-2))
    error('the mask file does not exist: %s', mask2(1:end-2))
end

% Check the number of Matlab workers open
if useparallel.do && isunix
    [~, w] = unix('top -n 1 | grep MATLAB');
    Nopen = numel(strfind(w, 'MATLAB'));
    [~, w] = unix('nproc');
    Ncpu = str2num(w);
    if Nopen > Ncpu
        fprintf('\n ***************************************************')
        fprintf('\n 	              WARNING                          ')
        fprintf('\n     There are more Matlab (%d) than cores (%d)     ', Nopen, Ncpu)
        fprintf('\n   The parallelization is likely to be inefficient  ')
        fprintf('\n ***************************************************')
        fprintf('\n')
    end
end
