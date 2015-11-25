% Script to create the subject-level PPI batch (Physiophysio):
% - take previously defined ROI (by loading the coordinates in the voxel
%   space)
% - extract the 1st eigenvalue of the BOLD signal in the ROIs with
%   spm_regions. This depends on a particular model, from which the effects
%   included in the design matrix (e.g. stimulus onsets, movement ...) are
%   regressed out, and whitening and filtering is re-used.
% - compute the PPI interaction term, on the deconvolved neural time series
% - specify the PPI specification, estimation and contrast batch. The GLM
%   include not only the interaction, main effects and session-specific
%   constant, but also the BOLD signal from seeds in the CSF and white
%   matter.
%
% can we m20150603_DistributePPI_computation.m to run the batch in parallel.

% =========================================================================
%                       REMOVE PREVIOUS VOI / PPI .mat
% =========================================================================
fprintf('\nCleaning...')
for iSub = 1:numel(sublist)
    
    % subject directory
    fname = sprintf('%s/subj%02.0f/first_level_estimates/%s', rootdir, sublist(iSub), modelname);
    
    % remove previous PPI and VOI file (if any)
    cd(fname)
    unix('rm -rf PPI*.mat');
    unix('rm -rf VOI*.mat');
end

fprintf('\n Extraction of signal in ROI...')
% loop over subject
for iSub = 1:numel(sublist)
    fprintf('\nsubject %d / %d', iSub, numel(sublist))
    
    % loop over ROI
    for iROI = 1:numel(ROIlist)+2;
        
        % =========================================================================
        %                               DEFINE ROI
        % =========================================================================

        fprintf('\nROI %d / %d', iROI, numel(ROIlist)+2)
        
        % Select ROI
        clear ClusterCoord_vx ROIname
        if iROI <= 2    % 'functional' ROI
            
            ClusterCoord_vx = ROIdef.ClusterCoord_vx{ROIlist(iROI)};
            ROIname         = ROIdef.ROIname{ROIlist(iROI)};
            
        else            % 'confound' ROI
            
            % Load subject and model specific SPM.mat of the estimated model
            % -> used to find the list of files (= scans)
            % -> used to find the precomputed filter and whitening matrix
            fname       = sprintf('%s/subj%02.0f/first_level_estimates/%s', rootdir, sublist(iSub), modelname);
            tmp         = load([fname, '/SPM.mat']);
            SPM         = tmp.SPM;
            
            % take the 1st EPI to have the actual vx / mni space of the data
            fname = deblank(SPM.xY.P(1,:));
            
            if iROI == 3 % take bilateral seed in the ventricules (3 mm radius)
                center_mm = [
                    -17 -35 18
                    17 -35 18];
                [~, ~, XYZvx_1] = GetData_SphereROI_nii(fname, center_mm(1,:), 3);
                [~, ~, XYZvx_2] = GetData_SphereROI_nii(fname, center_mm(2,:), 3);
                ClusterCoord_vx = [XYZvx_1, XYZvx_2];
                ROIname = 'Ventricule';
            end
            
            if iROI == 4 % take bilateral seed in the white matter (3 mm radius)
                center_mm = [
                    -28 -30 32
                    28 -30 32];
                [~, ~, XYZvx_1] = GetData_SphereROI_nii(fname, center_mm(1,:), 3);
                [~, ~, XYZvx_2] = GetData_SphereROI_nii(fname, center_mm(2,:), 3);
                ClusterCoord_vx = [XYZvx_1, XYZvx_2];
                ROIname = 'White Matter';
            end
            
            % check ROIs
            % tmp = spm_vol(fname);
            % fmspm12batch_CheckROI(ClusterCoord_vx, ROIname, tmp.mat)
        end
        
        % =========================================================================
        %                       EXTRACT SIGNAL FROM ROI
        % =========================================================================
        
        % Load subject and model specific SPM.mat of the estimated model
        % -> used to find the list of files (= scans)
        % -> used to find the precomputed filter and whitening matrix
        fname       = sprintf('%s/subj%02.0f/first_level_estimates/%s', rootdir, sublist(iSub), modelname);
        tmp         = load([fname, '/SPM.mat']);
        SPM         = tmp.SPM;
        
        % specific to this model: the model name 'Model01_QualityCheck' was renamed
        % into 'Model1_QualityCheck'. So, udpate SPM.swd
        SPM.swd = fname;
        
        % Move to the directory, so that the beta.nii can be read by spm_regions
        cd(fname)
        
        % Mask the ROI coordinates by the NaNs in the beta.nii (there are some at
        % near the cortex boundaries).
        flist       = dir('beta*.nii');
        betaval     = spm_get_data(flist(1).name, ClusterCoord_vx);
        ClusterCoord_vx = ClusterCoord_vx(:,~isnan(betaval));
        
        % Build the xY structure to specify ROI extraction options
        clear xY
        xY.Sess     = Inf;      % extract from all sessions
        xY.Ic       = NaN;      % NaN to remove all effects in the design matrix. Otherwise, specify the contrast indices (in SPM.xCon) to keep (specified as F-contrast), everything else is removed.
        xY.name     = sprintf('subj%02.0f_%s-nullspace_%s', sublist(iSub), modelname, ROIname); % name of the VOI
        xY.def      = 'all';    % take all voxels specified in the xSPM.XYZmm (hence, in ClusterCoord_vx)
        xY.xyz      = [];       % voxel coordinates in mm (left empty since all voxels specified are used
        
        % save VOIname
        VOIname{iSub, iROI} = xY.name;
        
        % Specify the VOI volume for get spm_regions in xSPM
        clear xSPM
        xSPM.XYZ    = ClusterCoord_vx; % used by spm_regions to extract data
        xSPM.XYZmm  = M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))]; % convert vx into mm (see spm_getSPM l.887)
        xSPM.M      = M;        % actually not used in spm_regions, but required.
        
        % Extract with spm_regions (automatically saved)
        % The 3rd argument, hReg, is left empty to avoid interaction with the GUI
        % (SPM graphical interface should be closed).
        [Y,xY] = spm_regions(xSPM, SPM, [], xY);
        
    end
end
fprintf('\n Extraction of signal in ROI completed')

% =========================================================================
%                       COMPUTE PPI
% =========================================================================
fprintf('\n\n Compute PPI deconvolution')

% loop over subject
for iSub = 1:numel(sublist)
    fprintf('\nsubject %d / %d', iSub, numel(sublist))
    
    % Specify model name
    % This is used to identify the output directory, and to get the recording
    % details (RT, ...)
    SPMname     = sprintf('%s/subj%02.0f/first_level_estimates/%s/SPM.mat', rootdir, sublist(iSub), modelname);
    
    cd(sprintf('%s/subj%02.0f/first_level_estimates/%s/', rootdir, sublist(iSub), modelname))
    
    % loop over sessions
    for iSess = 1:4
        
        % specify VOIs
        clear VOI
        VOI{1} = sprintf('VOI_%s_%d.mat', VOIname{iSub,1}, iSess);
        VOI{2} = sprintf('VOI_%s_%d.mat', VOIname{iSub,2}, iSess);
        VOI = char(VOI);
        
        % Specify output name (per session)
        ppiname = sprintf('subj%02.0f_PPI_Sess%d%', sublist(iSub), iSess);
        
        % Call the spm_peb_ppi routine
        % PPI = spm_peb_ppi(SPMname, ppiflag, VOI, Uu, ppiname, showGraphics)
        % Uu can be left empty since there is no psycho factor here
        % TO BE IMPROVED: the interaction term is absurd for me (low surp &
        % low conf are negative, thus result in a high value for the 
        % interaction term!!), fix it line 399 (post detrend: offset the
        % values so that there are all positive).
        PPI = spm_peb_ppi(SPMname, 'phipi', VOI, [], ppiname, 0);
        
        % NB: spm_peb_ppi filters and whiten the regressors using the same matrix
        % as for the 1st level estimation
        
    end
end

% CLEAN UP
% Move the VOI & PPI .mat file into a new directory.
for iSub = 1:numel(sublist)
    % create a 'PPI' folder
    target_dir = sprintf('%s/subj%02.0f/first_level_estimates/%s_%s/', ...
        rootdir, sublist(iSub), modelname, PPIname);
    mkdir(target_dir)
    
    orign_dir = sprintf('%s/subj%02.0f/first_level_estimates/%s/', ...
        rootdir, sublist(iSub), modelname);
    
    % Move VOI mat files
    flist = cellstr(spm_select('List', orign_dir, '^VOI.*\.mat'));
    if length(flist) == 1 && strcmp(flist{1}, '')
        % nothing to move
    else
        for iFile = 1:length(flist)
            movefile([orign_dir, flist{iFile}], [target_dir, flist{iFile}])
        end
    end
    
    % Move PPI mat files
    flist = cellstr(spm_select('List', orign_dir, '^PPI.*\.mat'));
    if length(flist) == 1 && strcmp(flist{1}, '')
        % nothing to move
    else
        for iFile = 1:length(flist)
            movefile([orign_dir, flist{iFile}], [target_dir, flist{iFile}])
        end
    end
end


% =========================================================================
%                       CREATE BATCH
% =========================================================================
fprintf('\n\n Write PPI batch')

% loop over subjects
for iSub = 1:numel(sublist)
    fprintf('\nsubject %d / %d', iSub, numel(sublist))
    
    % model name (fullpath of the folder for this subject)
    sub_model_name = sprintf('%s/subj%02.0f/first_level_estimates/%s_%s/', ...
        rootdir, sublist(iSub), modelname, PPIname);
    
    % Load subject and model specific SPM.mat of the estimated model
    % -> used to find the list of files (= scans)
    fname       = sprintf('%s/subj%02.0f/first_level_estimates/%s', rootdir, sublist(iSub), modelname);
    tmp         = load([fname, '/SPM.mat']);
    SPM         = tmp.SPM;
    
    % specify design matrix
    clear PPIdat
    desmat = zeros(sum(SPM.nscan), 3+2+numel(SPM.nscan)-1);
    scanindex = [0, cumsum(SPM.nscan)];
    for iSess = 1:4
        % load PPI data for this subject and this session
        PPIdat = load(sprintf('%s/PPI_subj%02.0f_PPI_Sess%d', ...
            sub_model_name, sublist(iSub), iSess));
        
        desmat((scanindex(iSess)+1):scanindex(iSess+1), 1)       = PPIdat.PPI.ppi;  % convolved interaction
        desmat((scanindex(iSess)+1):scanindex(iSess+1), 2)       = PPIdat.PPI.Y;    % BOLD region #1
        desmat((scanindex(iSess)+1):scanindex(iSess+1), 3)       = PPIdat.PPI.P;    % BOLD region #2
        
        % load confound (CSF and WS) for this subject and this session
        VOIdat = load(sprintf('%s/VOI_subj%02.0f_%s-nullspace_Ventricule_%d', ...
            sub_model_name, sublist(iSub), modelname, iSess));
        desmat((scanindex(iSess)+1):scanindex(iSess+1), 4)       = VOIdat.Y;        % CSF signal
        VOIdat = load(sprintf('%s/VOI_subj%02.0f_%s-nullspace_White Matter_%d', ...
            sub_model_name, sublist(iSub), modelname, iSess));
        desmat((scanindex(iSess)+1):scanindex(iSess+1), 5)       = VOIdat.Y;        % WM signal
        
        if iSess < numel(SPM.nscan)
            % a global constant is added by SPM, so having a constant for
            % the last session is useless
            desmat((scanindex(iSess)+1):scanindex(iSess+1), 5+iSess) = 1;           % constant
        end
    end
    R = desmat; % the variable must be called R to be read by SPM
    save([sub_model_name, 'desmat.mat'], 'R')
    
    % SPECIFY THE BATCH
    clear matlabbatch
    % Specification of the design matrix
    matlabbatch{1}.spm.stats.fmri_spec.dir              = {sub_model_name};
    
    matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'scans';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = SPM.xY.RT;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = 8;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans        = cellstr(SPM.xY.P);    % list of files
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond         = ...                   % leave empty
        struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi        = {''};                 % leave empty
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress      = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg    = {[sub_model_name, 'desmat.mat']};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf          = hpf;
    
    matlabbatch{1}.spm.stats.fmri_spec.fact              = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs  = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt              = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global            = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh           = spm_get_defaults('mask.thresh');
    matlabbatch{1}.spm.stats.fmri_spec.mask              = {mask};
    matlabbatch{1}.spm.stats.fmri_spec.cvi               = spm_get_defaults('stats.fmri.cvi');
    
    % Estimation the design matrix
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1)          = {sprintf('%s/SPM.mat', sub_model_name)};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals    = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical   = 1;                    % classical t-test
    
    % Report the contrast
    matlabbatch{3}.spm.stats.con.spmmat(1)               = {sprintf('%s/SPM.mat', sub_model_name)};
    
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = '(+) corr with interaction between seeds';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name    = '(-) corr with interaction between seeds';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 0 0 0 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name    = sprintf('(+) corr with %s', ROIdef.ROIname{ROIlist(1)});
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1 0 0 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name    = sprintf('(-) corr with %s', ROIdef.ROIname{ROIlist(1)});
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1 0 0 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name    = sprintf('(+) corr with %s', ROIdef.ROIname{ROIlist(2)});
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 1 0 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name    = sprintf('(-) corr with %s', ROIdef.ROIname{ROIlist(2)});
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 -1 0 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name    = sprintf('(+) corr with CSF');
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 1 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.name    = sprintf('(+) corr with WM');
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 1 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.delete                  = 1;
    
    
    % save batch
    save([sub_model_name, 'batch_PPI.mat'], 'matlabbatch')
end


