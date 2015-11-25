% Script to compare the residuals of two SPM models. 
% The paired-difference is computed at the subject level, the resulting 
% images are then smoothed and taken into a classical 2nd level, t-test
% analysis in SPM.
%
% Note that this is a quick 'model comparison'. As a BMS, it tests the
% goodness of the fit, not the consistency of the fitted parameter (in
% particular their direction) across subjects. Unlike BMS however, this
% script does not penalize for complexity, so that the comparison is fair 
% only if the number of regressors in the two design matrices are 
% identical.

% initialize the SPM & window.
clear all
try spm('ver');
catch
    error('the SPM toolbox is not in the Matlab''s path')
end
spm_results_ui('clear')
spm_results_ui('close')
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_path = spm('dir');

% GET DATA
% =========================================================================

% Get the rfx model directory (to get the model name)
models = cellstr(spm_select(2,'dir','Select 2 rfx models'));
modeldir1 = models{1};
modeldir2 = models{2};

% Process names
ind = strfind(modeldir1, '/');
modelname1 = modeldir1(ind(end)+1:end);
ind = strfind(modeldir2, '/');
modelname2 = modeldir2(ind(end)+1:end);
ind = strfind(modeldir1, '/group_analysis/');
rootdir = modeldir1(1:ind-1);

% Make short names
ind = strfind(modelname1, '_');
short_modelname1 = modelname1(1:ind(1)-1);
ind = strfind(modelname2, '_');
short_modelname2 = modelname2(1:ind(1)-1);


% Get list of subjects to plot
sublist = spm_input('subject list','+1','e','');
nSub = numel(sublist);

% Built list of individual ResMS.nii
ResMS1 = [];
ResMS2 = [];
for iSub = 1:nSub
    ResMS1(iSub,:) = ...
        sprintf('%s/subj%02.0f/first_level_estimates/%s/ResMS.nii,1', ...
        rootdir, sublist(iSub), modelname1);
    ResMS2(iSub,:) = ...
        sprintf('%s/subj%02.0f/first_level_estimates/%s/ResMS.nii,1', ...
        rootdir, sublist(iSub), modelname2);
end


% PERFORM THE 1ST LEVEL: PAIRED DIFFERENCE
% =========================================================================

% Use ImCalc to compute the paired difference of residual MS
for iSub = 1:nSub
    
    % Make the batch
    clear matlabbatch
    matlabbatch{1}.spm.util.imcalc.input = {char(ResMS1(iSub,:)) ; char(ResMS2(iSub,:))};
    matlabbatch{1}.spm.util.imcalc.output = sprintf('tmp_CompareRes_sub_%03.0f', iSub);
    matlabbatch{1}.spm.util.imcalc.outdir = {'/tmp/'};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1 - i2';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0; % NaNs are treated as zeros
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    % run the batch
    spm_jobman('run', matlabbatch);
end


% Get the output list of temporary files
List_diff = [];
for iSub = 1:nSub
    List_diff(iSub,:) = sprintf('/tmp/tmp_CompareRes_sub_%03.0f.nii,1', iSub);
end
List_diff = cellstr(char(List_diff));

% Smooth the images
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data      = { List_diff{:} }';
matlabbatch{1}.spm.spatial.smooth.fwhm      = [8 8 8];    % smoothing kernel
matlabbatch{1}.spm.spatial.smooth.dtype     = 0;                % output format = input
matlabbatch{1}.spm.spatial.smooth.im        = 0;                % no implicit masking
matlabbatch{1}.spm.spatial.smooth.prefix    = 's';              % new file prefix
spm_jobman('run', matlabbatch);


% PERFORM THE 2ND LEVEL: T-TEST
% =========================================================================

% Get smoothed files names
List_diff = cell(nSub,1);
for iSub = 1:nSub
    List_diff{iSub} = sprintf('/tmp/stmp_CompareRes_sub_%03.0f.nii,1', iSub);
end

% Parameters for the 2nd level analysis
dataloc.con_dir = sprintf('%s/group_analysis/Compare_ResMS/%s_vs_%s', ...
    rootdir, short_modelname1, short_modelname2); 
if ~exist(dataloc.con_dir, 'dir'), mkdir(dataloc.con_dir), end
param.mask = sprintf('%s/tpm/mask_ICV.nii,1', spm_path);
param.delprevcon = 1;
dataloc.con_name = sprintf('MSE %s > %s', short_modelname1, short_modelname2);


clear matlabbatch
[matfile, matlabbatch] = fmspm12batch_2ndlevel_WriteBatch_Classicalttest(dataloc, param, List_diff);
spm_jobman('run', matlabbatch);

% CLEAN UP
% =========================================================================
% remove the temporary files
! rm -f stmp_CompareRes_sub_*.nii
! rm -f tmp_CompareRes_sub_*.nii
