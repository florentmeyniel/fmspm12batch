% Script to set the experiment-specific pamareters for the pre-processing 
% of fMRI data.
% THIS IS AN EXAMPLE SCRIPT

% For the batch system to work properly, the folder architecture should be 
% as follows:
% /subjdir                  % name specified in this parameter file        
%   /funcdir                % name specified in this parameter file        
%   /anatdir                % name specified in this parameter file        
% Note that subjdir is datadir/subjXX; datadir is defined below as the data
% root directory and XX is the subject number, e.g. 1, 2, ... 20

% Specification of parameters
% =========================================================================

% Option for distributed computing
useparallel.do   = 0;               % 1 to run in parallel, 0 to run serially
useparallel.max  = 6;               % maximum number of batch launched in parallel
useparallel.cmd  = 'matlab-R2015a'; % command to invoke matlab from a terminal

% fsl command
flscmd           = 'fsl5.0';        % fls command use to call fsl functions, e.g. fsl5.0-fslroi

% subject list
sublist          = 1;
nSub             = length(sublist);

% actions to perform
% Possible actions are: 'standard' (slice timing, realign & unwrap based on
% a B0 file, i.e. a field map correction, normalize anat, corregister EPI
% and anat, normalize EPI, smooth EPI), 'AddTopup' (must be run after or
% together with 'standard': it stard from the realign & unwrap files and
% apply a second unwrapping for distortion due to the fast acquisition with
% AP/PA distortions, then corregister EPI with anat, normalize and smooth),
% 'run' to run the specified batch.
% Note that without 'run', the batch is simply saved.
actions          = {'run', 'SliceTiming', standard', 'Topup'};
% slicetiming or nothing
% realign or unrwrap
% topup or nothing
% seg/normalize
% smooth

% locate the data
spm_path         = '/volatile/meyniel/toolbox/matlab/spm12';
datadir          = '/neurospin/unicog/protocols/IRMf/Berkovitch_syntax_fMRI_2016/MRI_data/raw_data';
regexp_func      = '^epi.*\.nii';                    % regular expression to recognize functional sessions to analyze 
regexp_anat      = '^anat.*\.nii';                   % regular expression to recognize T1
funcdir          =  'fMRI';                          % path of fMRI data (4D nifti) within subject directory
anatdir          =  'anat';                          % path of anatomical image within subject directory
regexp_topupref  = '^uaepi_sess1_.*\.nii';           % regular expression to recognize the reference session used by FSL to realign the AP / PA volumes

% acquisition parameters
B0_TE            = [3.02 5.46];                      % short and long TE, in ms, of the B0 acquisition. (leave empty if no BO file)

% pre-processing parameters
smoothing_kernel = [5 5 5];                          % 1st level smoothing

% Automatically get info from the header
% =========================================================================
% Note on slice timing 
% --------------------
% SPM offers the possibility to specify directly the timing of each slice
% (which can be read from the dicom header). In this case, the acquisition
% time is not necessary (TA can be left to 0). This option is used in the
% current batch to deal with the multiband acquisition.
% Alternatively, on could specify the slice order (but this is not suitable
% for multiband acquisitions).
% Note by Cyril Poupon about the slice order on for interleaved acquisition
% on Siemens machines: 
% If the number of slice locations is even,
%   the slices are acquired as 2, 4, 6, ..., N, 1, 3, 5, ..., N-1.
% If the number of slice locations is odd, 
%   the slices are acquired as 1, 3, 5, ..., N-1, 2, 4, 6, ..., N. 

% get the slice timing info (specific to Neurospin).
% alternaltive: load them from a pre-existing file or set manually
% load timing info (generated by GetSliceTiming)
slice_timing = cell(nSub, 1);
TR = zeros(nSub, 1);
nslices = zeros(nSub, 1);
xyz_resol = zeros(nSub, 3);
currdir = pwd;
for iSub = 1:nSub
    % move to the subject's directory
    subjdir = sprintf('%s/subj%02.0f/%s/', datadir, iSub, funcdir);
    cd(subjdir)
    
    if ~exist('SliceTimingInfo.mat', 'file')
        % read information from the DICOM header on the server
        [slice_timing{iSub}, TR(iSub), ~, SliceThickness, ~, nslices(iSub), PixelSpacing] = ...
        fmspm12batch_preproc_GetSliceTiming_NS(sublist(iSub), funcdir, datadir, spm_path);
        xyz_resol(iSub) = [PixelSpacing(1), PixelSpacing(2), SliceThickness];   
    else
        % read the information previously saved
        tmp = load('SliceTimingInfo.mat');
        slice_timing{iSub} = tmp.SliceTiming;
        TR(iSub) = tmp.TR;
        nslices(iSub) = tmp.nslices;
        xyz_resol(iSub) = [tmp.PixelSpacing(1), tmp.PixelSpacing(2), tmp.SliceThickness];   
    end
end
cd(currdir)

% Check whether all participants have the same TR
if length(unique(TR)) == 1
    % all participants have the same TR => set it in a 1-value vector
    TR = unique(TR);
else
    error('different subjects have different TR')
end

% Check whether all participants have the same slice timing
if nSub > 1
    L = zeros(nSub, 1);
    for iSub = 1:nSub; L(iSub) = length(slice_timing{iSub}); end
    if length(unique(L)) == 1
        cat_ST = zeros(nSub, length(slice_timing{iSub}));
        for iSub =1:nSub; cat_ST(iSub,:) = slice_timing{iSub}; end
        diff_ST = diff(cat_ST);
        
        % check is the slice timing info are consistent at the 5ms resolution
        comparison = diff_ST > 5;
        if any(comparison(:))
            error(sprintf(['the acquisition parameters seem different over subjects...', ...
                '\n slice timing deviates >5ms for at least one slice!']))
        else
            % all participants the same slice timing => set it in a vector
            slice_timing = cat_ST(1,:);
        end
    else
        error(sprintf(['the acquisition parameters seem different over subjects...', ...
            '\n at least, the number of slide per volume differ!']))
    end
else
    slice_timing = slice_timing{1};
end

% check whether all participants have the same number of slices
if length(unique(nslices)) == 1
    % all participants have the same number of slices => set it in a 1-value vector
    nslices = unique(nslices);
else
    error('different subjects have different number of slices')
end

% check that resolution is the same for all subjects
if size(unique(xyz_resol, 'rows'), 1) == 1;
    xyz_resol = unique(xyz_resol, 'rows');
else
    error('different subjects have different resolutions')
end

% check that effective echo time are the same in all subjects
if (length(unique(total_readout_time_spm)) == 1) && (length(unique(total_readout_time_fsl)) == 1)
    total_readout_time_spm = unique(total_readout_time_spm);
    total_readout_time_fsl = unique(total_readout_time_fsl);
else
    error('different subjects have effective reaout time of EPIs')
end