function fmspm12batch_preproc(varargin)
% Function to write preprocessing batches, as specified by the input
% script that define parameters.
% This function also run (in parallel if specified) the batch.

% quit spm if open
try spm quit; end

% Get the parameters
if ~ischar(varargin{1})
    error('the argument should be a string')
else
    paramfile = varargin{1};
end
if ~strcmp(paramfile, 'RecursiveMode')
    % load parameters from file
    eval(sprintf('%s_paramfile', paramfile));
    
    dataloc.spm_path                = spm_path;            % spm directory
    dataloc.datadir                 = datadir;             % root directory for subject data
    dataloc.regexp_func             = regexp_func;         % regular expression to recognize functional sessions to analyze
    dataloc.regexp_anat             = regexp_anat;         % regular expression to recognize T1
    dataloc.funcdir                 = funcdir;             % path of fMRI data (4D nifti) within subject directory
    dataloc.anatdir                 = anatdir;             % path of anatomical image within subject directory
    
    param.TR                        = TR;
    param.nslices                   = nslices;             % number of slices
    param.voxel_size                = xyz_resol;
    param.smoothing_kernel          = smoothing_kernel;    % 1st level smoothing
    param.B0_TE                     = B0_TE;
    param.total_readout_time_spm    = total_readout_time_spm;
    param.total_readout_time_fsl    = total_readout_time_fsl;
    
    param.flscmd                    = flscmd;
    
    param.actions                   = lower(actions);
else
    % use parameter from input arguments
    param = varargin{2};
    dataloc = varargin{3};
end
sub_actions = {sub_actions{:}, 'slicetiming'}
sub_actions = [sub_actions, {'slicetiming'}]

% SPECIFICATION OF THE BATCH
% =========================================================================
if ismember('topup', actions)
    
    % DO ALL ACTIONS BEFORE TOPUP
    sub_actions = {};
    if ismember('slicetiming', actions); sub_actions = [sub_actions, {'slicetiming'}]; end
    if ismember('realign', actions);     sub_actions = [sub_actions, {'realign'}];     end
    if ismember('unwrap', actions);      sub_actions = [sub_actions, {'unwrap'}];      end
    if ~isempty(sub_actions)
        param.actions = {sub_actions{:}, 'run'};
        fmspm12batch_preproc('RecursiveMode', param, dataloc)
    end
    
    % DO THE TOPUP STEP
    prefix = '';
    if ismember('slicetiming', actions) && ismember('realign', actions)
        prefix = 'ra';
    end
    if ~ismember('slicetiming', actions) && ismember('realign', actions)
        prefix = 'r';
    end
    if ismember('slicetiming', actions) && ismember('unwrap', actions)
        prefix = 'ua';
    end
    if ~ismember('slicetiming', actions) && ismember('unwrap', actions)
        prefix = 'u';
    end
    dataloc.regexp_func = ['^', prefix, dataloc.regexp_func(2:end)];
    param.actions = {'AddTopupStep'};
    fmspm12batch_preproc('RecursiveMode', param, dataloc)
    
    % DO ALL ACTIONS AFTER TOPUP
    sub_actions = {};
    prefix = 't';
    if ismember('segmentnormalize', actions); sub_actions = [sub_actions, {'segmentnormalize'}]; end
    if ismember('smooth', actions);           sub_actions = [sub_actions, {'smooth'}];           end
    if ismember('slicetiming', actions) && ismember('realign', actions)
        prefix = 'tra';
    end
    if ~ismember('slicetiming', actions) && ismember('realign', actions)
        prefix = 'tr';
    end
    if ismember('slicetiming', actions) && ismember('unwrap', actions)
        prefix = 'tua';
    end
    if ~ismember('slicetiming', actions) && ismember('unwrap', actions)
        prefix = 'tu';
    end
    dataloc.regexp_func = ['^', prefix, dataloc.regexp_func(2:end)];
    if ~isempty(sub_actions)
        param.actions = {sub_actions{:}, 'run'};
        fmspm12batch_preproc('RecursiveMode', param, dataloc)
    end
    
else
    
    fprintf('\n Create batch for PREPROCESSING')    
    matfile = cell(nSub, 1); % full path of saved job
    job     = cell(nSub, 1); % job that is saved
    
    % create one job (matlabbatch) per subject
    for iSub = 1:nSub
        [matfile{iSub}, job{iSub}] = fmspm12batch_preproc_sf_make1job1sub(sublist(iSub), dataloc, param);
    end
    
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
end

% DEAL WITH EXTRA TOPUP CORRECTION
% =========================================================================
if ismember('AddTopupStep', actions)
    
    fprintf('\n Save Topup parameters for...')
    for iSub = 1:nSub
        fprintf('\n Subject %02.0f: %d/%d', sublist(iSub), iSub, nSub)
        
        fpath{iSub} = sprintf('%s/batch_specif_sub%02.0f_TopupCorrection.mat', ...
            pwd, sublist(iSub));
        SubNum = sublist(iSub);
        save(fpath{iSub}, 'SubNum', 'datadir', 'funcdir', ...
            'spm_path', 'regexp_topupref', 'total_readout_time_fsl');
    end
    
    if useparallel.do == 0
        
        % run batch serially
        fprintf('\n Run...')
        for iSub = 1:nSub
            
            fprintf('\n batch name: %s', fpath{iSub})
            
            % run batch
            fmspm12batch_AddTopupCorrection_job1sub(fpath{iSub})
        end
        
    elseif useparallel.do == 1
        
        % run batches in parallel
        fmspm12batch_runParalleljobs(fpath, dataloc.spm_path, ...
            useparallel.max, useparallel.cmd)
    end
end