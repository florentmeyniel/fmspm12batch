% Script to set the experiment-specific parameters for the pre-processing 
% of fMRI data.
% THIS IS AN EXAMPLE SCRIPT

% For the batch system to work properly, the folder architecture should be 
% as follows:
% /subjdir                  % name specified in this parameter file        
%   /funcdir                % name specified in this parameter file        
%   /anatdir                % name specified in this parameter file        
% Note that subjdir is datadir/subjXX; datadir is defined below as the data
% root directory and XX is the subject number, e.g. 01, 02, ... 20
%
% funcfir => where fMRI data can be found (.nii) with the prefix specified
% in this parameter file. 
% To model physiological artifacts, funcdir should also contains several
% files. There are 3 files per fMRI session: 
%   XXX_Info.log        % info from the scanner
%   XXX_PULS.log        % cardiac data
%   XXX_RESP.log        % respiratory data
% XXX is the name of the corresponding fMRI session: XXX.nii

% Specification of parameters
% =========================================================================

% Option for distributed computing
useparallel.do   = 0;               % 1 to run in parallel, 0 to run serially
useparallel.max  = 6;               % maximum number of batch launched in parallel
useparallel.cmd  = 'matlab';        % command to invoke matlab from a terminal

% fsl command
flscmd           = 'fsl5.0';        % fsl command to call fsl functions, e.g. fsl5.0-fslroi

% subject list
sublist          = 1;
nSub             = length(sublist);

% Actions to perform
% possible actions are: 
%   'retroicor'                                             no prefix (save text files)
%   'slicetiming'                                           prefix 'a'
%   'realign' or 'unrwarp' (the latter requires B0 files)   prefix 'r' and 'u' respectively
%   'topup' (requires AP/PA files)                          prefix 't'
%   'segmentnormalize'                                      prefix 'w'
%   'smooth'                                                prefix 's'
%   'run' to run the specified batch. Without 'run', the batch is simply saved (expect for topup)
% All these actions can be combined, e.g.: {'run', realign', 'segmentnormalize'}
% Because the batch works with any prefix defined by the user (see
% regexp_fun below), one can start the batch at any processing step. For
% instance, the slice timing correction may have been already computed, and
% the batch {'run', realign', 'segmentnormalize'} can start from there by
% specifying regexp_func = '^aepi.*\.nii';
% The order of arguments does not matter.
actions          = {'run', 'retroicor', 'slicetiming', 'realign', 'topup', 'segmentnormalize', 'smooth'};

% Locate the data
spm_path         = '/volatile/meyniel/toolbox/matlab/spm12';
datadir          = '/neurospin/unicog/protocols/IRMf/Berkovitch_syntax_fMRI_2016/MRI_data/raw_data/test_florent/';

regexp_func      = '^epi.*\.nii';                   % regular expression to recognize functional sessions to analyze (the ^ is mandatory)
regexp_anat      = '^anat.*\.nii';                  % regular expression to recognize T1 (the ^ is mandatory)
funcdir          =  'fMRI';                         % path of fMRI data (4D nifti) within subject directory
anatdir          =  'anat';                         % path of anatomical image within subject directory
regexp_topupref  = '^raepi_sess1_.*\.nii';          % functional session onto which field map files are aligned
topup_root_name  = 'ep2d';                          % the AP/PA files for topup are named 'topup_root_name_AP*.nii' 'topup_root_name_PA*.nii'

% acquisition parameters
B0_TE            = [];                              % short and long TE, in ms, of the B0 acquisition. (leave empty if no BO file)
blipdir          = 1;                               % phase encoding direction of EPI images (e.g. +1 if P->A; -1 is A->P)

% pre-processing parameters
smoothing_kernel = [5 5 5];                         % 1st level smoothing

% Topup options
topupOptions.do_realign      = 1;                   % realign AP/PA file onto reference EPI
topupOptions.do_estimate     = 1;                   % estimate the deformation field map
topupOptions.do_apply        = 1;                   % apply unwrapping
topupOptions.par.estimatemax = 6;                   % limit parallel process for estimate
topupOptions.par.applymax    = 3;                   % limit parallel process for apply (take lots of RAM)

% Retroicor options
retroicorOptions.toolbox     = '/volatile/meyniel/toolbox/matlab/PhysIO/'; % path to TAPAS toolbox

% Retroicor options for advanced users
retroicorOptions.verbose     = 0; % 0 = none; 1 = main plots (default);  2 = debugging plots, for setting up new study; 3 = all
retroicorOptions.order.c     = 3; % order of the Fourier expansion for the cardiac phase
retroicorOptions.order.r     = 4; % order of the Fourier expansion for the respiratory phase
retroicorOptions.order.cr    = 1; % order of the Fourier expansion for the interaction cardiac x respiratory 

% Automatically get acquisition parameters from the DICOM header
% =========================================================================
% Note on slice timing 
% --------------------
% SPM offers the possibility to specify directly the timing of each slice
% (which can be read from the DICOM header). In this case, the acquisition
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

% fmspm12batch_preproc_GetSliceTiming_NS reads the acquisition parameters
% from the header of the DICOM. By default, the functions search for the
% relevant DICOM files on the Neurospin server, but the full path name of
% the file can also be provided directly.
slice_timing            = cell(nSub, 1);
TR                      = zeros(nSub, 1);
nslices                 = zeros(nSub, 1);
xyz_resol               = zeros(nSub, 3);
total_readout_time_spm  = zeros(nSub, 1);
total_readout_time_fsl  = zeros(nSub, 1);
currdir = pwd;
for iSub = 1:nSub
    % move to the subject's directory
    subjdir = sprintf('%s/subj%02.0f/%s/', datadir, sublist(iSub), funcdir);
    cd(subjdir)
    
    if ~exist('SliceTimingInfo.mat', 'file')
        % read information from the DICOM header on the server
        [slice_timing{iSub}, TR(iSub), ~, SliceThickness, ~, nslices(iSub), ...
            PixelSpacing, total_readout_time_spm(iSub), total_readout_time_fsl(iSub)] = ...
        fmspm12batch_preproc_GetSliceTiming_NS(sublist(iSub), funcdir, datadir, spm_path);
        xyz_resol(iSub,:) = [PixelSpacing(1), PixelSpacing(2), SliceThickness];   
    else
        % read the information previously saved
        tmp = load('SliceTimingInfo.mat');
        slice_timing{iSub} = tmp.SliceTiming;
        TR(iSub) = tmp.TR;
        nslices(iSub) = tmp.NumberOfSlices;
        xyz_resol(iSub,:) = [tmp.PixelSpacing(1), tmp.PixelSpacing(2), tmp.SliceThickness];  
        total_readout_time_spm(iSub) = tmp.total_readout_time_spm;
        total_readout_time_fsl(iSub) = tmp.total_readout_time_fsl;
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

% get that actions are recognized
actions = lower(actions);
for iAction = 1:numel(actions)
    if ~ismember(actions{iAction}, {'run', 'retroicor', 'slicetiming', 'realign', 'unwarp', 'topup', 'segmentnormalize', 'smooth'})
        error('action ''%s'' not recognized', actions{iAction})
    end
end
