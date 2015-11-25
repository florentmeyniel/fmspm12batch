function fmspm12batch_2ndlevel(modelname)
% Function to write the second-level batches, as specified by the input
% script that define parameters.
% -> can smooth the 1st level contrast
% -> can perform a 2nd level t-test.
%
% syntax: fmspm12batch_2ndlevel('mymodel') will callmymodel_paramfile.m to get
% the parameters and know what to do.


% INITIALIZATION
% =========================================================================

tic
% quit spm if open
try spm quit; end

% Get the parameters
eval(sprintf('%s_paramfile', modelname));

dataloc.spm_path       = spm_path;            % spm directory
dataloc.datadir        = datadir;             % root directory for subject data
dataloc.modelname      = modelname;           % get model name
dataloc.funcdir        = funcdir ;            % subfolder were EPIs are

param.mask             = mask2;               % mask for estimation
param.skernek          = smoothing_kernel2;   % smoothing kernel for 2nd level analysis
param.delprevcon       = deleteprevcon;       % delete previous contrast
param.nSub             = nSub;
param.contrast_type    = contrast_type;	      % con or scon
param.sublist          = sublist;


% Make directories to save results
if ~exist(sprintf('%s/group_analysis/rfxmodel/%s/', datadir, modelname), 'dir')
    mkdir(sprintf('%s/group_analysis/rfxmodel/%s/', datadir, modelname))
end


% SMOOTHING FOR SECOND LEVEL ANALYSIS
% =========================================================================
if ismember('smooth', actions)
    
    fprintf('\n Smooth 1st level contrast images: [%d %d %d] mm', param.skernek)
    
    % Get list of contrasts, if not specified
    if isempty(con_list)
        % get SPM.mat of 1st subject
        SPM = load(sprintf('%s/subj%02.0f/first_level_estimates/%s/SPM.mat', ...
                datadir, sublist(1), modelname));
            
        % get number of contrast in this subject (assuming it is a good
        % reference for the other subjects)
        con_list = 1:numel(SPM.SPM.xCon);
    end
    
    % List contrast to smooth
    n = 0;
    for iSub = 1:param.nSub
        for iCon = 1:numel(con_list)
            n = n + 1;
            confilelist{n} = ...
                sprintf('%s/subj%02.0f/first_level_estimates/%s/con_%04.0f.nii,1', ...
                dataloc.datadir, param.sublist(iSub), dataloc.modelname, con_list(iCon));
        end
    end
    
    matlabbatch{1}.spm.spatial.smooth.data      = confilelist';     
    matlabbatch{1}.spm.spatial.smooth.fwhm      = param.skernek;    % smoothing kernel
    matlabbatch{1}.spm.spatial.smooth.dtype     = 0;                % output format = input
    matlabbatch{1}.spm.spatial.smooth.im        = 0;                % no implicit masking
    matlabbatch{1}.spm.spatial.smooth.prefix    = 's';              % new file prefix
    
    % Save the SPM batch
    matfile = sprintf('%s/group_analysis/rfxmodel/%s/batch_smooth.mat', datadir, modelname);
    if exist(matfile)
        delete(matfile); fprintf('\nremove previous batch\n')
    end
    
    fprintf('\n Writing SPM batch at: \n %s\n', matfile)
    save(matfile,'matlabbatch');
end

if ismember('smooth', actions) && ismember('run', actions)
    fmspm12batch_run1job(matfile, dataloc.spm_path)
end


% SECOND LEVEL T-TEST
% =========================================================================

if ismember('rfx_ttest', actions)
    
    % Print the list of subjects
    fprintf('\n List of subjects taken in the 2nd level:')
    for iSub = 1:nSub
        fprintf('\n   subj%02.0f', sublist(iSub))
    end
    fprintf('\n Contrast type: %s', contrast_type)
    fprintf('\n')
    
    % Get list of contrasts, if not specified
    if isempty(con_list)
        % get SPM.mat of 1st subject
        SPM = load(sprintf('%s/subj%02.0f/first_level_estimates/%s/SPM.mat', ...
                datadir, sublist(1), modelname));
            
        % get number of contrast in this subject (assuming it is a good
        % reference for the other subjects)
        con_list = 1:numel(SPM.SPM.xCon);
    end
 
    % Initialize jobs
    nCon    = numel(con_list);
    matfile = cell(nCon, 1); % full path of saved job
    job     = cell(nCon, 1); % job that is saved
    
    for iCon = 1:nCon
        
        fprintf('\n Contrast %d: ', iCon)
        
        % Make some checks
        % ~~~~~~~~~~~~~~~~
        con_name = cell(1, nSub);
        for iSub = 1:nSub
            
            % Check that the contrast exists
            con_fullpath = sprintf('%s/subj%02.0f/first_level_estimates/%s/%s_%04.0f.nii', ...
                datadir, sublist(iSub), modelname, contrast_type, con_list(iCon));
            if ~exist(con_fullpath, 'file')
                error('This constrast does not exist: \n %s', con_fullpath)
            end
            
            % Get contrast names
            SPM = load(sprintf('%s/subj%02.0f/first_level_estimates/%s/SPM.mat', ...
                datadir, sublist(iSub), modelname));
            
            con_name{iSub} = SPM.SPM.xCon(con_list(iCon)).name;
        end
        
        % Check that contrast names are the same for all subjects
        for iSub = 2:nSub
            if ~strcmp(con_name{1}, con_name{iSub})
                error('(s)con_%04.0f has different names: \n subj%02.0f: %s \n subj%02.0f: %s', ...
                    con_list(iCon), sublist(1), con_name{1}, sublist(iSub), con_name{iSub})
            end
        end
        
        % Move on with this contrast name
        tmp = con_name{1}; clear con_name; con_name = tmp; clear tmp;
        fprintf('%s', con_name)
        
        % Write the 2nd level batch
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Make the directory to save the results for this contrast
        con_dir = sprintf('%s/group_analysis/rfxmodel/%s/%s', datadir, modelname, con_name);
        if exist(con_dir, 'dir')
            if deleteprevcon == 1
                unix(sprintf('rm -rf ''%s''', con_dir))
                mkdir(con_dir)
            end
        else
            mkdir(con_dir)
        end
        
        % Make list of files for the contrast
        dataloc.con_dir = con_dir;
        dataloc.con_name = con_name;
        confilelist = cell(param.nSub,1);
        for iSub = 1:param.nSub
            confilelist{iSub} = sprintf('%s/subj%02.0f/first_level_estimates/%s/%s_%04.0f.nii,1', ...
                dataloc.datadir, param.sublist(iSub), dataloc.modelname, param.contrast_type, con_list(iCon));
        end
                
        % write the SPM batch
        [matfile{iCon}, job{iCon}] = fmspm12batch_2ndlevel_WriteBatch_Classicalttest(dataloc, param, confilelist);

    end
end

if ismember('rfx_ttest', actions) && ismember('run', actions)
    
    fprintf('\n Run batch for CONTRAST 2ND LEVEL T-TEST')
    
    if useparallel.do == 0
        
        % run batch serially
        for iCon = 1:length(con_list)
            fprintf('\n batch name: %s', matfile{iCon})
            fmspm12batch_run1job(matfile{iCon}, dataloc.spm_path)
        end
        
    elseif useparallel.do == 1
        
        % run batches in parallel
        fmspm12batch_runParalleljobs(matfile, dataloc.spm_path, ...
            useparallel.max, useparallel.cmd)
    end
end

fprintf('\n End of the 2nd level analysis\n')
toc

