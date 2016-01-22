function [matfile, matlabbatch] = fmspm12batch_1stlevel_specify_sf_make1job1sub(iSub, dataloc, param)
% Function to specify the 1st-level specification of the design matrix.
%
% [matfile, matlabbatch] = fmspm12batch_1stlevelspecify_sf_make1job1sub(iSub, dataloc, param)
%  	Outputs:
% 		* matfile: full path of the batch saved.
% 		* matlabbatch: batch structure (that is saved)
% 	Inputs:
% 		* iSub: subject number, to find the subject folder subjXX
% 		* dataloc: structure to locate data
% 		* param: structure with specification parameters.
%
% This function will look for .mat where regressors are specified.
% These .mat should be found at:
%  %s/subj%d/MutliCond/%s_multicond_session%d.mat; dataloc.datadir, iSub, dataloc.modelname, session #

% unpack input data
% =========================================================================
spm_path         = dataloc.spm_path;            % spm directory  
datadir          = dataloc.datadir;             % root directory for subject data: dataloc/subjXX
regexp_func      = dataloc.regexp_func;         % regular expression to recognize functional sessions to analyze 
funcdir          = dataloc.funcdir;             % path of fMRI data (4D nifti) within subject directory
modelname        = dataloc.modelname;           % name of the model
TR               = param.TR;                    % TR in s
bases_functions  = param.bases_functions;       % temporal basis functions
refslice         = param.refslice;              % reference slice for timing within a volume
mask             = param.mask;                  % mask for estimation

% Format the temporal basis function in the SPM nomenclature
switch bases_functions
    case 'hrf' % HRF simple
        basisfunc = [0 0];
    case 'hrf+deriv'
        basisfunc = [1 0];
    case 'hrf+2derivs'
        basisfunc = [1 1];
    case 'fir' % FIR
        basisfunc = struct('length',14.4,'order',12);
    case 'fourier'
        basisfunc = struct('length',14.4,'order',12);
    otherwise
        error('Unknown basis function');
end


% initialize spm
% =========================================================================
addpath(spm_path)
spm('defaults', 'FMRI');


% Get the data files
% =========================================================================

% Get subject directory (with anat & fMRI)
subjdir = sprintf('%s/subj%02.0f/', datadir, iSub);

% Get functional files to analyze
fdir =  sprintf('%s%s/', subjdir, funcdir);
ffiles = spm_select('List', fdir, regexp_func);
nrun = size(ffiles,1);
if nrun == 0
    warning(sprintf('No functional file found for %s', ...
        subjects{csubj}))
    return
end
funcfiles = cell(nrun, 1);
cffiles = cellstr(ffiles);
for i = 1:nrun
    funcfiles{i} = cellstr(spm_select('ExtFPList', fdir, ['^', cffiles{i}], Inf));
end

% get where to save the SPM.mat (full path)
modeldir = sprintf('%sfirst_level_estimates/%s', subjdir, modelname);

% GET THE MOVEMENT PARAMETER FILES
% check that epi_ was used as a prefix during the importation of the data
for i = 1:nrun
    if isempty(strfind(cffiles{i}, 'epi_'))
        error('cannot find the prefix epi_ in the file name %s', cffiles{i})
    end
end

% get all realignemnt parameters
all_movpar = cellstr(spm_select('FPList', fdir, '^rp_.*\.txt'));

% match file name
movpar = cell(nrun,1);
for i = 1:nrun
    % get the unique file identifier
    start_ind = strfind(cffiles{i}, 'epi_') + 4;
    FileUniqueTag = cffiles{i}(start_ind:end-4);
    
    % match the identifier between the realigement and functional files
    MatchedFileInd = cellfun(@(x) ~isempty(strfind(x, FileUniqueTag)), all_movpar);
    movpar{i} = all_movpar{MatchedFileInd};
end

% Specify design matrix
% =========================================================================
matlabbatch{1}.spm.stats.fmri_spec.dir 				= {modeldir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units 		= 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT 			= TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t 		= spm_get_defaults('stats.fmri.t');
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 		= refslice; % reference time bin for microtime resolution.
for iRun = 1:nrun
    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).scans 	= funcfiles{iRun};
    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond 		= struct(...
									'name', {}, ...
									'onset', {}, ...
									'duration', {}, ...
									'tmod', {}, ...
									'pmod', {}, ...
									'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi 	= {sprintf('%sMultiCond/%s_multicond_session%d.mat', ...
									subjdir, modelname, iRun)};
    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).regress 	= struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi_reg 	= {movpar{iRun}};
    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).hpf 		= spm_get_defaults('stats.fmri.hpf');
end
matlabbatch{1}.spm.stats.fmri_spec.fact 			= struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = basisfunc;
matlabbatch{1}.spm.stats.fmri_spec.volt 			= 1;
matlabbatch{1}.spm.stats.fmri_spec.global 			= 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh 			= spm_get_defaults('mask.thresh');
matlabbatch{1}.spm.stats.fmri_spec.mask 			= {mask};
matlabbatch{1}.spm.stats.fmri_spec.cvi 				= spm_get_defaults('stats.fmri.cvi');


% Save SPM batch
% =========================================================================
matfile = sprintf('%s/batch_1stlevel_specify.mat', modeldir);
if exist(matfile) 
    delete(matfile); fprintf('\nremove previous batch\n')
end
if ~exist(modeldir); mkdir(modeldir); end

fprintf('\n Writing SPM batch for subject %d:', iSub)
fprintf('\n %s\n', matfile)
save(matfile,'matlabbatch');

