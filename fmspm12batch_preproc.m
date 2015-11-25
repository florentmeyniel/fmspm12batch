function fmspm12batch_preproc(paramfile)
% Function to write preprocessing batches, as specified by the input 
% script that define parameters.
% This function also run (in parallel if specified) the batch.

% quit spm if open
try spm quit; end

% Get the parameters
eval(paramfile);

dataloc.spm_path       = spm_path;            % spm directory  
dataloc.datadir        = datadir;             % root directory for subject data
dataloc.regexp_func    = regexp_func;         % regular expression to recognize functional sessions to analyze 
dataloc.regexp_anat    = regexp_anat;         % regular expression to recognize T1
dataloc.funcdir        = funcdir;             % path of fMRI data (4D nifti) within subject directory
dataloc.anatdir        = anatdir;             % path of anatomical image within subject directory

param.nslices          = nslices;               % number of slices
param.deltaEPI         = deltaEPI;              % readout between 2 EPI in ms ('Ecart echo' in the Siemens PDF)
param.iPAT             = iPAT;                  % EPI acceleration
param.voxel_size       = voxel_size; 
param.smoothing_kernel = smoothing_kernel;      % 1st level smoothing

% SPECIFICATION OF THE BATCH
% =========================================================================

matfile = cell(nSub, 1); % full path of saved job
job     = cell(nSub, 1); % job that is saved
% create one job (matlabbatch) per subject
fprintf('\n Create batch for PREPROCESSING')
for iSub = 1:nSub
    if length(TR) > 1
        param.TR = TR(iSub);
    else
        param.TR = TR;
    end
    fprintf('\n subject %d: TR=%d s', sublist(iSub), param.TR)
    if isstruct(slice_timing)
        param.slice_timing = slice_timing{iSub};
    else
        param.slice_timing = slice_timing;
    end
    
    [matfile{iSub}, job{iSub}] = fmspm12batch_preproc_sf_make1job1sub(sublist(iSub), dataloc, param);
end

% EXECUTION OF THE BATCH
% =========================================================================

if ismember('run', actions)
    fprintf('\n Run batch for PREPROCESSING')
    
    if useparallel.do == 0
        
        % run batch serially
        for iSub = 1:nSub
            
            fprintf('\n batch name: %s', matfile{iSub})
            
            % run batch
            fmspm12batch_run1job(matfile{iSub}, dataloc.spm_path)
        end
        
    elseif useparallel.do == 1
        
        % run batches in parallel
        fmspm12batch_runParalleljobs(matfile, dataloc.spm_path, ...
            useparallel.max, useparallel.cmd)
    end
end
