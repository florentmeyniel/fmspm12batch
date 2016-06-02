function fmspm12batch_AddTopupCorrection_job1sub(speciffile)
% Function to add a topup correction just after the 'Realign and Unwrap'
% images computed in the fmspm12batch_preproc pipeline. The corrected
% images have the 't' in their prefix.
% The function also run the normalization (wtua*) and smoothing (swtua*),
% as in the standard pipeline.
% NB: it should take 1h30 per subject.

% load specification file
sp = load(speciffile);

SubNum                  = sp.SubNum;
datadir                 = sp.datadir;
funcdir                 = sp.funcdir;
spm_path                = sp.spm_path;
regexp_func             = sp.regexp_func;
regexp_topupref         = sp.regexp_topupref;
total_readout_time_fsl  = sp.total_readout_time_fsl;

% Get the specific info on this subject
% initialization
% =========================================================================
% initialize SPM
addpath(spm_path)
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Get subject directory (with anat & fMRI)
subjdir = sprintf('%s/subj%02.0f/', datadir, SubNum);
fdir =  sprintf('%s%s/', subjdir, funcdir);

% Identify the AP PA file for calibrating the deformation, with the nii
% extention ,1
B0_AP = spm_select('ExtFPList', fdir, '^ep2d_AP_.*\.nii', 1);
B0_PA = spm_select('ExtFPList', fdir, '^ep2d_PA_.*\.nii', 1);

% Identify the 1st session, after the unwrapping step
EPIref = spm_select('ExtFPList', fdir, regexp_topupref, 1);

% check that files exist
if strcmp(B0_AP, ''); error('cannot find B0_AP file'); end
if strcmp(B0_PA, ''); error('cannot find B0_PA file'); end
if strcmp(EPIref, ''); error('cannot find EPIref file'); end

% REALIGN THE AP & PA CALIBRATION FILE ONTO 1st EPI
% =========================================================================

% Make the SPM batch
clear estwrite matlabbatch
estwrite.data               = {{EPIref}, {B0_AP}, {B0_PA}}';
estwrite.eoptions.quality   = spm_get_defaults('realign.estimate.quality');                              
estwrite.eoptions.sep       = spm_get_defaults('realign.estimate.sep');
estwrite.eoptions.fwhm      = spm_get_defaults('realign.estimate.fwhm');
estwrite.eoptions.rtm       = 0;            % register to 1st
estwrite.eoptions.interp    = spm_get_defaults('realign.estimate.interp');
estwrite.eoptions.wrap      = [0 0 0];      % no wraping (because phase encoding is different!)
estwrite.eoptions.weight    = '';           % no weighting
estwrite.roptions.which     = [1 0];        % Do not reslice the 1st image
estwrite.roptions.interp    = spm_get_defaults('realign.write.interp');
estwrite.roptions.wrap      = [0 0 0];      % no wraping (because phase encoding is different!)
estwrite.roptions.mask      = 1;            % mask with voxel that are not 0s
estwrite.roptions.prefix    = 'r';          % prefix

matlabbatch{1}.spm.spatial.realign.estwrite = estwrite; 

% Run the job
spm_jobman('run', matlabbatch);

% ESTIMATE THE DEFORMATION MAP WITH FSL & TOPUP
% =========================================================================

% Merge the realigned AP & PA nii files in the same file with .nii.gz format
B0_AP = spm_select('List', fdir, '^rep2d_AP_.*\.nii');
B0_PA = spm_select('List', fdir, '^rep2d_PA_.*\.nii');
cmd = sprintf('cd %s; fsl5.0-fslmerge -t b0_APPA %s %s', ...
                fdir, B0_AP(1:end-4), B0_PA(1:end-4)); 
unix(cmd)

% Create a text file with the direction of phase encoding. 
% A>P is -1; P>A is 1 
% (could be checked with Romain's script topup_param_from_dicom). 
cmd = sprintf('cd %s; echo $''0 -1 0 %6.5f \n0 1 0 %6.5f'' > acq_param.txt', ...
                fdir, total_readout_time_fsl, total_readout_time_fsl);
unix(cmd)

% Compute deformation with Topup (this takes ~20 minutes and only 1 CPU)
% The result that will be used is APPA_DefMap
% For sanity checks, sanitycheck_DefMap is the deformation field and
% sanitycheck_unwarped_B0 are the corrected images
% NB: in the b02b0.cnf, subsamp starts =2, which implies that the size
% of the image cannot be an odd number. Here we force subsampling =1.
fprintf('\n Compute the APPA deformation with Topup')
cmd = sprintf(['cd %s; fsl5.0-topup ', ...
    '--imain=b0_APPA --datain=acq_param.txt --config=b02b0.cnf ', ...
    '--out=APPA_DefMap --fout=sanitycheck_DefMap --iout=sanitycheck_unwarped_B0 ', ...
    '--subsamp=1,1,1,1,1,1,1,1,1'], ...
    fdir);
unix(cmd)

% APPLY THE CORRECTION
% =========================================================================
fprintf('\n Apply the topup correction...')

% List the unwrapped EPI, with ua prefix
EPI = cellstr(spm_select('List', fdir, regexp_func));

for iFile = 1:numel(EPI)
    fprintf('\n  Sess: %d/%d', iFile, numel(EPI))
    % Make the topup command to apply the correction
    % the 't' prefix is added to the output file
    cmd = sprintf(['cd %s; fsl5.0-applytopup --imain=%s ', ...
        '--inindex=2 ', ...                                      % because the EPIs are in PA, the 2nd row of the acq_param.txt file
        '--topup=APPA_DefMap --datain=acq_param.txt ', ...
        '--out=%s --method=jac'], ...
        fdir, EPI{iFile}(1:end-4), ['t', EPI{iFile}(1:end-4)]);
    unix(cmd)
    
    % Conver the results from .nii.gz to .nii
    cmd = sprintf('cd %s; fsl5.0-fslchfiletype NIFTI %s', ...
        fdir, ['t', EPI{iFile}(1:end-4)]);
    unix(cmd)
end
fprintf('\n')
