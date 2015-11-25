function fmspm12batch_AddTopupCorrection(paramfile)
% Batch to add the Topup correction for the AP/PA deformation on the EPI
% images. The function loads the parameters from the preprocessing 
% parameter file, and loop over subject to apply the correction.
% This batch is like a 'patch' to the standard preprocessing: it should be
% run only after fmspm12batch_preproc.m


% quit spm if open
try spm quit; end

% Get the parameters
eval(paramfile);


for iSub = 1:nSub
    fprintf('\n Subject %02.0f: %d/%d', sublist(iSub), iSub, nSub)
    fmspm12batch_AddTopupCorrection_job1sub(...
        sublist(iSub), datadir, funcdir, anatdir, ...
        voxel_size, smoothing_kernel, spm_path)
end