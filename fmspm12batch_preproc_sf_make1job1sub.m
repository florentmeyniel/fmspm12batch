function [matfile, matlabbatch] = fmspm12batch_preproc_sf_make1job1sub(iSub, dataloc, param)

% Create and save an SPM batch, for one subject, with the following
% preprocessing steps:
% - Slice Timing
% - Estimation of FieldMap deformations (NB: mag & phase files should have specific names...)
% - Segment anatomy and align on an MNI template
% - Realign & Unwrap, wrt. first EPI
% - Coregister mean EPI on anatoy
% - Normalize the anatomy, given the normalized segmentation
% - Normalize the EPIs, given the normalized segmentation
% - Smooth EPIs

% =========================================================================
%                       UNPACK INPUT DATA
% =========================================================================
spm_path                = dataloc.spm_path;             % spm directory
datadir                 = dataloc.datadir;              % root directory for subject data
regexp_func             = dataloc.regexp_func;          % regular expression to recognize functional sessions to analyze
regexp_anat             = dataloc.regexp_anat;          % regular expression to recognize T1
funcdir                 = dataloc.funcdir;              % path of fMRI data (4D nifti) within subject directory
anatdir                 = dataloc.anatdir;              % path of anatomical image within subject directory

nslices                 = param.nslices;                % number of slices
voxel_size              = param.voxel_size;
smoothing_kernel        = param.smoothing_kernel;       % 1st level smoothing
TR                      = param.TR;
slice_timing            = param.slice_timing;
B0_TE                   = param.B0_TE;                  % short and long TE, in ms, of the B0 acquisition.
flscmd                  = param.flscmd;                 % fls command use to call fsl function, e.g. fsl5.0-fslroi
total_readout_time_spm  = param.total_readout_time_spm; % effective readout time of EPIs
total_readout_time_fsl  = param.total_readout_time_fsl; % effective readout time of EPIs

actions                 = param.actions;

% =========================================================================
%                           INITIALIZE SPM
% =========================================================================
addpath(spm_path)
spm('defaults', 'FMRI');
spm_jobman('initcfg');
stage = 0;

% =========================================================================
%                            GET DATA FILES
% =========================================================================

% Get subject directory (with anat & fMRI)
subjdir = sprintf('%s/subj%02.0f/', datadir, iSub);

% Get T1 file to analyze
adir =  sprintf('%s/%s/', subjdir, anatdir);
anatfile = spm_select('FPList', adir, regexp_anat);
if isequal(anatfile,  '')
    warning(sprintf('No anat file found for %s', iSub))
    return
end

% Get functional files to analyze
fdir =  sprintf('%s%s/', subjdir, funcdir);
ffiles = spm_select('List', fdir, regexp_func);
nrun = size(ffiles,1);
if nrun == 0
    warning(sprintf('No functional file found for %s', iSub))
    return
end
funcfiles = cell(nrun, 1);
cffiles = cellstr(ffiles);
for i = 1:nrun
    funcfiles{i} = cellstr(spm_select('ExtFPList', fdir, ['^', cffiles{i}], Inf));
end

% =========================================================================
%                              SLICE TIMING
% =========================================================================
if ismember('slicetiming', actions)
    stage = stage + 1;
    
    matlabbatch{stage}.spm.temporal.st.scans = { funcfiles{:} };
    matlabbatch{stage}.spm.temporal.st.nslices = nslices;
    matlabbatch{stage}.spm.temporal.st.tr = TR;
    matlabbatch{stage}.spm.temporal.st.ta = 0;              % TA can be left to 0 if slice timing is explit (as ms, not order)
    matlabbatch{stage}.spm.temporal.st.so = slice_timing;
    matlabbatch{stage}.spm.temporal.st.refslice = 0;        % reference slice is actually ref. time (in ms) because slice times are provided.
    matlabbatch{stage}.spm.temporal.st.prefix = 'a';
end

% =========================================================================
%                   RELALIGN AND UNWRAP WITH A FIELD MAP
% =========================================================================
if ismember('unwrap', actions)
    stage = stage + 1;
    stage_fieldmap = stage;
    
    % Compute FieldMap correction
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % arrange B0 files and identify magnitude and phase (should be checked
    % manually at the importation of the data!!)
    if exist([fdir, 'B0'], 'dir');
        unix(sprintf('rm -rf %s', [fdir, 'B0']))
        mkdir([fdir, 'B0'])
    else
        mkdir([fdir, 'B0'])
    end
    fname_B0_1 = spm_select('List', fdir, '^B0_1.*\.nii');
    copyfile([fdir, fname_B0_1], strcat([fdir, 'B0/', 'B0_mag_', fname_B0_1(6:end)]))
    
    fname_B0_2 = spm_select('List', fdir, '^B0_2.*\.nii');
    copyfile([fdir, fname_B0_2], strcat([fdir, 'B0/', 'B0_phase_', fname_B0_2(6:end)]))
    
    % extract the short & intermediate TE with FSL
    % 1st (0) = shortest: 3.00 ms
    fname_B0mag = spm_select('List', [fdir, 'B0'], '^B0_mag.*\.nii');
    cmd = sprintf('cd %s; %s-fslroi %s %s 0 1', ...
        [fdir, 'B0'], flscmd, ...
        fname_B0mag, strcat([fname_B0mag(1:7), 'shortTE_', fname_B0mag(8:end)]));
    unix(cmd)
    fname_B0phase = spm_select('List', [fdir, 'B0'], '^B0_phase.*\.nii');
    cmd = sprintf('cd %s; %s-fslroi %s %s 0 1', ...
        [fdir, 'B0'], flscmd, ...
        fname_B0phase, strcat([fname_B0phase(1:9), 'shortTE_', fname_B0phase(10:end)]));
    unix(cmd)
    
    % 2nd (1) = intermediate: 5.46 ms
    cmd = sprintf('cd %s; %s-fslroi %s %s 1 1', ...
        [fdir, 'B0'], flscmd, ...
        fname_B0mag, strcat([fname_B0mag(1:7), 'longTE_', fname_B0mag(8:end)]));
    unix(cmd)
    cmd = sprintf('cd %s; %s-fslroi %s %s 1 1', ...
        [fdir, 'B0'], flscmd, ...
        fname_B0phase, strcat([fname_B0phase(1:9), 'longTE_', fname_B0phase(10:end)]));
    unix(cmd)
    
    % convert these files to nii
    unix(sprintf('cd %s; fsl5.0-fslchfiletype NIFTI B0_mag_shortTE_%s'  , ...
        [fdir, 'B0'], fname_B0mag(8:end)))
    unix(sprintf('cd %s; fsl5.0-fslchfiletype NIFTI B0_mag_longTE_%s'   , ...
        [fdir, 'B0'], fname_B0mag(8:end)))
    unix(sprintf('cd %s; fsl5.0-fslchfiletype NIFTI B0_phase_shortTE_%s', ...
        [fdir, 'B0'], fname_B0phase(10:end)))
    unix(sprintf('cd %s; fsl5.0-fslchfiletype NIFTI B0_phase_longTE_%s' , ...
        [fdir, 'B0'], fname_B0phase(10:end)))
    
    % get names of files
    fname_shortmag   = spm_select('FPList', [fdir, '/B0'], '^B0_mag_shortTE.*\.nii');
    fname_longmag    = spm_select('FPList', [fdir, '/B0'], '^B0_mag_longTE.*\.nii');
    fname_shortphase = spm_select('FPList', [fdir, '/B0'], '^B0_phase_shortTE.*\.nii');
    fname_longphase  = spm_select('FPList', [fdir, '/B0'], '^B0_phase_longTE.*\.nii');
    
    % Create batch to compute the Voxel Deformation Map (VDM)
    % Note that one VDM is computed per session (and aligned on the 1st EPI of
    % the session)
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.shortphase = {fname_shortphase};
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.shortmag   = {fname_shortmag};
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.longphase  = {fname_longphase};
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.longmag    = {fname_longmag};
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.et = B0_TE;                        % [shortTE longTE] in ms
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.maskbrain = 1;
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.blipdir = 1;                       % +1 for P -> A (would be -1 for A -> P)
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.tert = total_readout_time_spm;     % EPI readout time
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.epifm = 0;                         % 0: fieldmap is not an EPI image (it is a gradiant echo)
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.ajm = 0;                           % no jacobian modulation (as recommanded by SPM)
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.method = pm_get_defaults('UNWRAPPING_METHOD');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.fwhm = pm_get_defaults('FWHM');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.pad = pm_get_defaults('PAS');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.ws = pm_get_defaults('WS');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.template = pm_get_defaults('MFLAGS.TEMPLATE');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.fwhm = pm_get_defaults('MFLAGS.FWHM');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.nerode = pm_get_defaults('MFLAGS.NERODE');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.ndilate = pm_get_defaults('MFLAGS.NDILATE');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.thresh = pm_get_defaults('MFLAGS.THRESH');
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.reg = pm_get_defaults('MFLAGS.REG');
    for irun = 1:nrun
        matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.session(irun).epi = {funcfiles{irun}{1}};           % specify 1st EPI for session alignement (for quality control with a display)
    end
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.matchvdm = 1;                                           % so that the VDM is aligned on the 1st EPI of each session
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.sessname = 'session';                                   % VDM file extension
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.writeunwarped = 1;
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.anat = {[anatfile, ',1']};
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.matchanat = 1;                                          % for visual comparison
    
    
    % Realign and unwrap the functional images
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Note that in SPM, if multiple sessions are provided, the 1st scan of each
    % session are realigned to each other (on the first scan of the first session)
    % and then, within each session, the scan are realigned to the 1st scan of that
    % session. The realignment parameters saved are thus relative to the 1st scan
    % of the 1st session, which is not a problem since the 1st scan of the different
    % sessions are first aligned to each other.
    stage = stage + 1;
    
    if ismember('slicetiming', actions)
        prefix = 'a';
    else
        prefix = '';
    end
    for iRun = 1:nrun
        % add prefix to file, depending on previous preprocessing steps
        afflist = AddPrefix(fflist, prefix);
        
        % list files for each session
        matlabbatch{stage}.spm.spatial.realignunwarp.data(iRun).scans = afflist{iRun};
        matlabbatch{stage}.spm.spatial.realignunwarp.data(iRun).pmscan(1) = ...
            cfg_dep(sprintf('Phase and Magnitude Data: Voxel displacement map (Subj 1, Session %d)', iRun), ...
            substruct('.','val', '{}',{stage_fieldmap}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('()',{1}, '.','vdmfile', '{}',{iRun}));
    end
    
    matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.quality = spm_get_defaults('realign.estimate.quality');
    matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.sep = spm_get_defaults('realign.estimate.sep');
    matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.fwhm = spm_get_defaults('realign.estimate.fwhm');
    matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.rtm = spm_get_defaults('realign.estimate.rtm');
    matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.einterp = spm_get_defaults('realign.estimate.interp');
    matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.ewrap = [0 1 0];                                      % wraping along the Y axis (because acquisition in PA direction)
    matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.weight = '';
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.basfcn = spm_get_defaults('unwarp.estimate.basfcn');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.regorder = spm_get_defaults('unwarp.estimate.regorder');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.lambda = spm_get_defaults('unwarp.estimate.regwgt');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.jm = spm_get_defaults('unwarp.estimate.jm');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.fot = spm_get_defaults('unwarp.estimate.foe');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.sot = spm_get_defaults('unwarp.estimate.soe');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.uwfwhm = spm_get_defaults('unwarp.estimate.fwhm');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.rem = spm_get_defaults('unwarp.estimate.rem');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.noi = spm_get_defaults('unwarp.estimate.noi');
    matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.expround = spm_get_defaults('unwarp.estimate.expround');
    matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.uwwhich = spm_get_defaults('realign.write.which');
    matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.rinterp = spm_get_defaults('realign.write.interp');
    matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.wrap = [0 1 0]; % wraping along the Y axis (because acquisition in PA direction)
    matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.mask = spm_get_defaults('realign.write.mask');
    matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
end

% Simple re-alignment
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ismember('realign', actions)
    
    stage = stage + 1;
    
    % get prefix of files, depending on previous processing steps
    if ismember('slicetiming', actions)
        prefix = '^a';
    else
        prefix = '^';
    end
    
    for iRun = 1:nrun
        % add prefix to file, depending on previous preprocessing steps
        afflist = AddPrefix(funcfiles, prefix);
        
        % list files for each session
        matlabbatch{stage}.spm.spatial.realign.estwrite.data{iRun} = afflist{iRun};
    end
    
    
    
    matlabbatch{stage}.spm.spatial.realign.estwrite.eoptions.quality = spm_get_defaults('realign.estimate.quality');
    matlabbatch{stage}.spm.spatial.realign.estwrite.eoptions.sep = spm_get_defaults('realign.estimate.quality');
    matlabbatch{stage}.spm.spatial.realign.estwrite.eoptions.fwhm = spm_get_defaults('realign.estimate.fwhm');
    matlabbatch{stage}.spm.spatial.realign.estwrite.eoptions.rtm = spm_get_defaults('realign.estimate.rtm');               % register to 1st
    matlabbatch{stage}.spm.spatial.realign.estwrite.eoptions.interp = spm_get_defaults('realign.estimate.interp');
    matlabbatch{stage}.spm.spatial.realign.estwrite.eoptions.wrap = [0 1 0];        % wraping along the Y axis (because acquisition in PA direction)
    matlabbatch{stage}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{stage}.spm.spatial.realign.estwrite.roptions.which = spm_get_defaults('realign.write.which');
    matlabbatch{stage}.spm.spatial.realign.estwrite.roptions.interp = spm_get_defaults('realign.write.interp');
    matlabbatch{stage}.spm.spatial.realign.estwrite.roptions.wrap = [0 1 0];        % wraping along the Y axis (because acquisition in PA direction)
    matlabbatch{stage}.spm.spatial.realign.estwrite.roptions.mask = spm_get_defaults('realign.write.mask');
    matlabbatch{stage}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
end

% =========================================================================
%                   SEGMENT AND NORMALIZE THE DATA
% =========================================================================


if ismember('segmentnormalize', actions)
    
    % Segmentation and spatially normalize anatomy with an MNI template by
    % computing an affine transformation
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % (using SPM default parameters)
    
    stage = stage + 1;
    stage_segmentation = stage;
    
    matlabbatch{stage}.spm.spatial.preproc.channel.vols = { anatfile };
    matlabbatch{stage}.spm.spatial.preproc.warp.write = [0 1];
    matlabbatch{stage}.spm.spatial.preproc.channel.biasreg = 0.001; % light regularisation (suggest by default)
    matlabbatch{stage}.spm.spatial.preproc.channel.biasfwhm = 60; % cutoff in mm (suggest by default)
    matlabbatch{stage}.spm.spatial.preproc.channel.write = [0 1]; % save bias corrected image
    ngaus  = [1 1 2 3 4 2];
    native = [1 1 1 0 0 0];
    for c = 1:6 % tissue class c (values suggest by default)
        matlabbatch{stage}.spm.spatial.preproc.tissue(c).tpm = {
            fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
        matlabbatch{stage}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
        matlabbatch{stage}.spm.spatial.preproc.tissue(c).native = [native(c) 0];
        matlabbatch{stage}.spm.spatial.preproc.tissue(c).warped = [0 0];
    end
    matlabbatch{stage}.spm.spatial.preproc.warp.mrf = 1; % (suggest by default)
    matlabbatch{stage}.spm.spatial.preproc.warp.cleanup = 1; % (suggest by default)
    matlabbatch{stage}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2]; % (suggest by default)
    matlabbatch{stage}.spm.spatial.preproc.warp.affreg = 'mni'; % (suggest by default)
    matlabbatch{stage}.spm.spatial.preproc.warp.fwhm = 0; % (suggest by default)
    matlabbatch{stage}.spm.spatial.preproc.warp.samp = 3; % (suggest by default)
    matlabbatch{stage}.spm.spatial.preproc.warp.write = [1 1]; % write the forward & inverse deformation fields
    
    % Get directory of bias corrected image to use for brain image
    stage = stage + 1;
    stage_braindir = stage;
    matlabbatch{stage_braindir}.cfg_basicio.file_dir.cfg_fileparts.files(1) = ...
        cfg_dep('Segment: Bias Corrected (1)', ...
        substruct('.','val', '{}',{stage_segmentation}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    
    
    % Create a brain image (skull-stripped bias corrected) with Image Calculator
    stage = stage + 1;
    stage_Imcal = stage;
    for c = 1:3
        matlabbatch{stage_Imcal}.spm.util.imcalc.input(c) = ...
            cfg_dep(sprintf('Segment: c%d Images', c), ...
            substruct('.','val', '{}',{stage_segmentation}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{c}, '.','c', '()',{':'}));
    end
    matlabbatch{stage_Imcal}.spm.util.imcalc.input(4) = ...
        cfg_dep('Segment: Bias Corrected (1)', ...
        substruct('.','val', '{}',{stage_segmentation}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{stage_Imcal}.spm.util.imcalc.output = 'Brain';
    matlabbatch{stage_Imcal}.spm.util.imcalc.outdir(1) = ...
        cfg_dep('Get Pathnames: Directories (unique)', ...
        substruct('.','val', '{}',{stage_braindir}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','up'));
    matlabbatch{stage_Imcal}.spm.util.imcalc.expression = '(i1 + i2 + i3) .* i4';
    matlabbatch{stage_Imcal}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{stage_Imcal}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{stage_Imcal}.spm.util.imcalc.options.mask = 0;
    matlabbatch{stage_Imcal}.spm.util.imcalc.options.interp = 1;
    matlabbatch{stage_Imcal}.spm.util.imcalc.options.dtype = 4;
    
    
    % coregistration: 1st EPI to Brain
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % (using SPM default parameters)
    
    % get the calculated anatomical brain image for this subject
    ImCalBrain = [adir, 'Brain.nii,1'];
    
    % get EPIs
    % List the EPI, with ua prefix
    if  ismember('slicetiming', actions) && ismember('realign', actions) &&  ismember('topup', actions); prefix = 'tra'; end
    if ~ismember('slicetiming', actions) && ismember('realign', actions) &&  ismember('topup', actions); prefix = 'tr';  end
    if  ismember('slicetiming', actions) && ismember('unwrap', actions)  &&  ismember('topup', actions); prefix = 'tua'; end
    if ~ismember('slicetiming', actions) && ismember('unwrap', actions)  &&  ismember('topup', actions); prefix = 'tu';  end
    if  ismember('slicetiming', actions) && ismember('realign', actions) && ~ismember('topup', actions); prefix = 'ra';  end
    if ~ismember('slicetiming', actions) && ismember('realign', actions) && ~ismember('topup', actions); prefix = 'r';   end
    if  ismember('slicetiming', actions) && ismember('unwrap', actions)  && ~ismember('topup', actions); prefix = 'ua';  end
    if ~ismember('slicetiming', actions) && ismember('unwrap', actions)  && ~ismember('topup', actions); prefix = 'u';   end
    OtherEPI = AddPrefix(funcfiles, prefix);
    
    stage = stage + 1;
    
    matlabbatch{stage}.spm.spatial.coreg.estimate.ref(1)            = {ImCalBrain};
    matlabbatch{stage}.spm.spatial.coreg.estimate.source(1)         = {EPIref};
    matlabbatch{stage}.spm.spatial.coreg.estimate.other             = OtherEPI;
    matlabbatch{stage}.spm.spatial.coreg.estimate.eoptions.cost_fun = spm_get_defaults('coreg.estimate.cost_fun');
    matlabbatch{stage}.spm.spatial.coreg.estimate.eoptions.sep      = spm_get_defaults('coreg.estimate.sep');
    matlabbatch{stage}.spm.spatial.coreg.estimate.eoptions.tol      = spm_get_defaults('coreg.estimate.tol');
    matlabbatch{stage}.spm.spatial.coreg.estimate.eoptions.fwhm     = spm_get_defaults('coreg.estimate.fwhm');
    
    
    % Spatial normalisation of EPIs, given the normalized segmentation (and its
    % affine transformation)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % (using SPM default parameters)
    
    stage = stage + 1;
    
    matlabbatch{stage}.spm.spatial.normalise.write.subj.def(1) = ... 
        cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{stage_segmentation}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
    matlabbatch{stage}.spm.spatial.normalise.write.subj.resample    = OtherEPI;
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb      = spm_get_defaults('normalise.write.bb');
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.vox     = voxel_size;
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.interp  = spm_get_defaults('normalise.write.interp');
    
    
    % Spatial normalization of anat, given the normalized segmentation (and its
    % affine transformation)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % (using SPM default parameters)
    
    stage = stage + 1;
    
    matlabbatch{stage}.spm.spatial.normalise.write.subj.def         = {ForwardDef};
    matlabbatch{stage}.spm.spatial.normalise.write.subj.resample(1) = {ImCalBrain};
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb      = spm_get_defaults('normalise.write.bb');
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.vox     = [1 1 1];
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.interp  = spm_get_defaults('normalise.write.interp');
end

% =========================================================================
%                               SMOOTH EPIs
% =========================================================================
if ismember('smooth', actions)
    stage = stage + 1;
    
    % get normalized EPIs
    if  ismember('slicetiming', actions) && ismember('realign', actions) &&  ismember('topup', actions); prefix = 'wtra'; end
    if ~ismember('slicetiming', actions) && ismember('realign', actions) &&  ismember('topup', actions); prefix = 'wtr';  end
    if  ismember('slicetiming', actions) && ismember('unwrap', actions)  &&  ismember('topup', actions); prefix = 'wtua'; end
    if ~ismember('slicetiming', actions) && ismember('unwrap', actions)  &&  ismember('topup', actions); prefix = 'wtu';  end
    if  ismember('slicetiming', actions) && ismember('realign', actions) && ~ismember('topup', actions); prefix = 'wra';  end
    if ~ismember('slicetiming', actions) && ismember('realign', actions) && ~ismember('topup', actions); prefix = 'wr';   end
    if  ismember('slicetiming', actions) && ismember('unwrap', actions)  && ~ismember('topup', actions); prefix = 'wua';  end
    if ~ismember('slicetiming', actions) && ismember('unwrap', actions)  && ~ismember('topup', actions); prefix = 'wu';   end
    NormEPIs = AddPrefix(funcfiles, prefix);
    
    matlabbatch{stage}.spm.spatial.smooth.data                      = NormEPIs;
    matlabbatch{stage}.spm.spatial.smooth.fwhm                      = smoothing_kernel;
    matlabbatch{stage}.spm.spatial.smooth.dtype                     = 0; % output data type = same as input
    matlabbatch{stage}.spm.spatial.smooth.im                        = 0; % no implicit masking
    matlabbatch{stage}.spm.spatial.smooth.prefix                    = 's';
end

% =========================================================================
%                            SAVE SPM BATCH
% =========================================================================
matfile = sprintf('%s/batch_preproc_sub%02.0f.mat', datadir, iSub);
if exist(matfile)
    delete(matfile); fprintf('\nremove previous batch\n')
end
fprintf('\n Writing SPM batch for subject %d:', iSub)
fprintf('\n %s\n', matfile)
save(matfile,'matlabbatch');

end
% ~~~~~~~~~~~~~~~~~~~~~~~  SUBFUNCTIONS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function afflist = AddPrefix(fflist, prefix)
% add a prefix to the file names in fflist

% initialize
afflist = cell(size(fflist));

for iiRun = 1:numel(fflist)
    
    if size(fflist{iiRun},1) > 1
        afflist{iiRun} = cell(size(fflist{iiRun}));
        for iiFile = 1:numel(fflist{iiRun})
            % get start of the file name
            ind = strfind(fflist{iiRun}{iiFile}, '/');
            
            if isempty(ind)
                % add prefix
                afflist{iiRun}{iiFile} = [prefix, fflist{iiRun}{iiFile}];
            else
                ind = max(ind);
                
                % add prefix
                afflist{iiRun}{iiFile} = ...
                    [fflist{iiRun}{iiFile}(1:ind), prefix, fflist{iiRun}{iiFile}(ind+1:end)];
            end
        end
    else
        % get start of the file name
        ind = strfind(fflist{iiRun}, '/');
        if isempty(ind)
            % add prefix
            afflist{iiRun} = [prefix, fflist{iiRun}];
        else
            ind = max(ind);
            % add prefix
            afflist{iiRun} = ...
                [fflist{iiRun}(1:ind), prefix, fflist{iiRun}(ind+1:end)];
        end
    end
end
end
