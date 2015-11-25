function fmspm12batch_1stlevel(modelname)
% Function to write the first-level batches, as specified by the input 
% script that define parameters.
%
% e.g. fmspm12batch_1stlevel('mymodel') will call:
%   - mymodel_paramfile.m to get the parameters and what to do
%   - mymodel_multicond.m to get the regressors
%   - mymodel_contrast.m to get the contrasts

% INITIALIZATION
% =========================================================================

tic
% quit spm if open
try spm quit; end

% Get the parameters
eval(sprintf('%s_paramfile', modelname));
addpath(spm_path)

dataloc.spm_path       = spm_path;            % spm directory  
dataloc.datadir        = datadir;             % root directory for subject data
dataloc.regexp_func    = regexp_func;         % regular expression to recognize functional sessions to analyze 
dataloc.modelname      = modelname;           % get model name
dataloc.funcdir        = funcdir ;            % subfolder were EPIs are

param.bases_functions  = bases_functions;     % temporal basis functions
param.refslice 	       = refslice;            % reference slice in a volume
param.mask             = mask;                % mask for estimation

% SPECIFICATION OF THE DESIGN MATRIX
% =========================================================================

% create one job (a matlabbatch.mat) per subject
if ismember('specify', actions)
    
    fprintf('\n Create batch for SPECIFICATION OF DESIGN MATRIX')
    
    % Get the 'multicond', which specifies the regressors
    eval(sprintf('%s_multicond', modelname));
    
    matfile = cell(nSub, 1); % full path of saved job
    job     = cell(nSub, 1); % job that is saved
    for iSub = 1:nSub
        
        fprintf('\n ... Subject %d, (ID: subj %02.0f)', iSub, sublist(iSub))

	% Get the acquisition TR
        if length(TR) > 1
            param.TR = TR(iSub);
        else
            param.TR = TR;
        end
        
        [matfile{iSub}, job{iSub}] = ...
            fmspm12batch_1stlevel_specify_sf_make1job1sub(sublist(iSub), dataloc, param);
    end
end

if ismember('specify', actions) && ismember('run', actions)
    
    fprintf('\n Run batch for SPECIFICATION OF DESIGN MATRIX')
    % NB: the specication calls spm_run_fmri_spec.m & spm_fMRI_design
    
    if useparallel.do == 0
        
        % run batch serially
        for iSub = 1:nSub
            
            fprintf('\n batch name: %s', matfile{iSub})
            
            % remove the SPM.mat (if any)
            dat  = load(matfile{iSub});
            fspm = sprintf('%s/SPM.mat', dat.matlabbatch{1}.spm.stats.fmri_spec.dir{1});
            if exist(fspm)
                fprintf('\n Remove previous SPM.mat')
                delete(fspm)
            end
            
            % run batch
            fmspm12batch_run1job(matfile{iSub}, dataloc.spm_path)
        end
        
    elseif useparallel.do == 1
        
        % run batches in parallel
        fmspm12batch_runParalleljobs(matfile, dataloc.spm_path, ...
            useparallel.max, useparallel.cmd)
    end
end

% ESTIMATION OF THE DESIGN MATRIX
% =========================================================================

if ismember('estimate', actions)
    
    fprintf('\n Create batch for ESTIMATION OF DESIGN MATRIX')
    
    matfile = cell(nSub, 1); % full path of saved job
    for iSub = 1:nSub
        
        % Get location of the SPM.mat specified
        SPMloc = sprintf('%s/subj%02.0f/first_level_estimates/%s/SPM.mat', ...
            datadir, sublist(iSub), modelname);
        
        if ~exist(SPMloc)
            error(['Could not find the specification of the 1st level:\n %s', ...
                '\n maybe the specification level has not been run...'], SPMloc)
        end
        
        % Specify batch
        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {SPMloc};    % full path of SPM.mat
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;  % save residuals
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % 1 for classical, 2 & 3 for Bayesian
        
        % save batch
        matfile{iSub} = ...
            sprintf('%s/subj%02.0f/first_level_estimates/%s/batch_1stlevel_estimate.mat', ...
            datadir, sublist(iSub), modelname);
        if exist(matfile{iSub})
            delete(matfile{iSub}); fprintf('\nremove previous batch\n')
        end
        
        fprintf('\n Writing SPM batch for estimation, 1st level, subject %d (ID: subj %02.0f)', iSub, sublist(iSub))
        fprintf('\n %s\n', matfile{iSub})
        save(matfile{iSub},'matlabbatch');
    end
end

if ismember('estimate', actions) && ismember('run', actions)

    fprintf('\n Run batch for ESTIMATION OF DESIGN MATRIX')
    
    if useparallel.do == 0
        
        % run batch serially
        for iSub = 1:nSub
            fprintf('\n batch name: %s', matfile{iSub})
            fmspm12batch_run1job(matfile{iSub}, dataloc.spm_path)
        end
        
    elseif useparallel.do == 1
        
        % run batches in parallel
        fmspm12batch_runParalleljobs(matfile, dataloc.spm_path, ...
            useparallel.max, useparallel.cmd)
    end
end

% SPECIFICATION AND ESTIMATION OF THE CONTRASTS
% =========================================================================
if ismember('contrast', actions)
    
    fprintf('\n Create batch for CONTRASTS')
    
    % Get the contrasts
    eval(sprintf('%s_contrast', modelname));
    
    matfile = cell(nSub, 1); % full path of saved job        
    for iSub = 1:nSub
        
        % Get location of the SPM.mat previously specified
        SPMloc = sprintf('%s/subj%02.0f/first_level_estimates/%s/SPM.mat', ...
            datadir, sublist(iSub), modelname);
        
        if ~exist(SPMloc)
            error(['Could not find the specification of the 1st level:\n %s', ...
                '\n maybe the specification level has not been run...'], SPMloc)
        end

        % Specify batch
        clear matlabbatch
        matlabbatch{1}.spm.stats.con.spmmat = {SPMloc};
        consess = fmspm12batch_1stlevel_contrast(datadir, sublist(iSub), modelname);
        matlabbatch{1}.spm.stats.con.consess = consess; 
        matlabbatch{1}.spm.stats.con.delete = deleteprevcon;

        % save batch
        matfile{iSub} = ...
            sprintf('%s/subj%02.0f/first_level_estimates/%s/batch_1stlevel_contrast.mat', ...
            datadir, sublist(iSub), modelname);
        if exist(matfile{iSub})
            delete(matfile{iSub}); fprintf('\nremove previous batch\n')
        end
        
        fprintf('\n Writing SPM batch for contrast, 1st level, subject %d (ID: subj %02.0f)', iSub, sublist(iSub))
        fprintf('\n %s\n', matfile{iSub})
        save(matfile{iSub},'matlabbatch');
    end
end

if ismember('contrast', actions) && ismember('run', actions)

    fprintf('\n Run batch for CONTRASTS')
    
    if useparallel.do == 0
        
        % run batch serially
        for iSub = 1:nSub
            fprintf('\n batch name: %s', matfile{iSub})
            fmspm12batch_run1job(matfile{iSub}, dataloc.spm_path)
        end
        
    elseif useparallel.do == 1
        
        % run batches in parallel
        fmspm12batch_runParalleljobs(matfile, dataloc.spm_path, ...
            useparallel.max, useparallel.cmd)
    end
end

fprintf('\n End of the 1st level analysis\n')
toc

