function [matfile, matlabbatch] = fmspm12batch_preproc_sf_make1job1sub(iSub, dataloc, param)

% Create and save an SPM batch, for one subject, with the following
% preprocessing steps:
% - Slice Timing
% - Estimation of FieldMap deformations (NB: mag & phase files should have specific names...)
% - Realign and unwrap, or simple realign, wrt. first EPI
% - Segment anatomy and align on an MNI template
% - Coregister mean EPI on anatomy
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
adir =  sprintf('%s%s/', subjdir, anatdir);
anatfile = spm_select('FPList', adir, regexp_anat);
if isequal(anatfile,  '')
    error('No anat file found for %s', iSub)
end

% Get functional files to analyze
fdir =  sprintf('%s%s/', subjdir, funcdir);
ffiles = spm_select('List', fdir, regexp_func);
nrun = size(ffiles,1);
if nrun == 0
    error('No functional file found for %s', iSub)
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
    clear st
    
    st.scans    = { funcfiles{:} };
    st.nslices  = nslices;
    st.tr       = TR;
    st.ta       = 0;                % TA can be left to 0 if slice timing is explit (as ms, not order)
    st.so       = slice_timing;
    st.refslice = 0;                % reference slice is actually ref. time (in ms) because slice times are provided.
    st.prefix   = 'a';
    
    matlabbatch{stage}.spm.temporal.st = st;
end

% =========================================================================
%                           REALIGN
% =========================================================================

% OPTION 1: realign and unrwrap with a field map
if ismember('unwarp', actions)
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
    clear subj
    subj.shortphase                             = {fname_shortphase};
    subj.shortmag                               = {fname_shortmag};
    subj.longphase                              = {fname_longphase};
    subj.longmag                                = {fname_longmag};
    subj.defaults.defaultsval.et                = B0_TE;                    % [shortTE longTE] in ms
    subj.defaults.defaultsval.maskbrain         = 1;
    subj.defaults.defaultsval.blipdir           = 1;                        % +1 for P -> A (would be -1 for A -> P)
    subj.defaults.defaultsval.tert              = total_readout_time_spm;   % EPI readout time
    subj.defaults.defaultsval.epifm             = 0;                        % 0: fieldmap is not an EPI image (it is a gradiant echo)
    subj.defaults.defaultsval.ajm               = 0;                        % no jacobian modulation (as recommanded by SPM)
    subj.defaults.defaultsval.uflags.method     = pm_get_defaults('UNWRAPPING_METHOD');
    subj.defaults.defaultsval.uflags.fwhm       = pm_get_defaults('FWHM');
    subj.defaults.defaultsval.uflags.pad        = pm_get_defaults('PAS');
    subj.defaults.defaultsval.uflags.ws         = pm_get_defaults('WS');
    subj.defaults.defaultsval.mflags.template   = pm_get_defaults('MFLAGS.TEMPLATE');
    subj.defaults.defaultsval.mflags.fwhm       = pm_get_defaults('MFLAGS.FWHM');
    subj.defaults.defaultsval.mflags.nerode     = pm_get_defaults('MFLAGS.NERODE');
    subj.defaults.defaultsval.mflags.ndilate    = pm_get_defaults('MFLAGS.NDILATE');
    subj.defaults.defaultsval.mflags.thresh     = pm_get_defaults('MFLAGS.THRESH');
    subj.defaults.defaultsval.mflags.reg        = pm_get_defaults('MFLAGS.REG');
    for irun = 1:nrun
        subj.session(irun).epi = {funcfiles{irun}{1}};                      % specify 1st EPI for session alignement (for quality control with a display)
    end
    subj.matchvdm = 1;                                                      % so that the VDM is aligned on the 1st EPI of each session
    subj.sessname = 'session';                                              % VDM file extension
    subj.writeunwarped = 1;
    subj.anat = {[anatfile, ',1']};
    subj.matchanat = 1;                                                     % for visual comparison
    
    matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj = subj;
    
    
    % Realign and unwrap the functional images
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Note that in SPM, if multiple sessions are provided, the 1st scan of each
    % session are realigned to each other (on the first scan of the first session)
    % and then, within each session, the scan are realigned to the 1st scan of that
    % session. The realignment parameters saved are thus relative to the 1st scan
    % of the 1st session, which is not a problem since the 1st scan of the different
    % sessions are first aligned to each other.
    stage = stage + 1;
    clear realignunwarp
    
    if ismember('slicetiming', actions)
        prefix = 'a';
    else
        prefix = '';
    end
    for iRun = 1:nrun
        % add prefix to file, depending on previous preprocessing steps
        afflist = AddPrefix(fflist, prefix);
        
        % list files for each session
        realignunwarp.data(iRun).scans = afflist{iRun};
        realignunwarp.data(iRun).pmscan(1) = ...
            cfg_dep(sprintf('Phase and Magnitude Data: Voxel displacement map (Subj 1, Session %d)', iRun), ...
            substruct('.','val', '{}',{stage_fieldmap}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('()',{1}, '.','vdmfile', '{}',{iRun}));
    end
    
    realignunwarp.eoptions.quality      = spm_get_defaults('realign.estimate.quality');
    realignunwarp.eoptions.sep          = spm_get_defaults('realign.estimate.sep');
    realignunwarp.eoptions.fwhm         = spm_get_defaults('realign.estimate.fwhm');
    realignunwarp.eoptions.rtm          = spm_get_defaults('realign.estimate.rtm');
    realignunwarp.eoptions.einterp      = spm_get_defaults('realign.estimate.interp');
    realignunwarp.eoptions.ewrap        = [0 1 0];                                      % wraping along the Y axis (because acquisition in PA direction)
    realignunwarp.eoptions.weight       = '';
    realignunwarp.uweoptions.basfcn     = spm_get_defaults('unwarp.estimate.basfcn');
    realignunwarp.uweoptions.regorder   = spm_get_defaults('unwarp.estimate.regorder');
    realignunwarp.uweoptions.lambda     = spm_get_defaults('unwarp.estimate.regwgt');
    realignunwarp.uweoptions.jm         = spm_get_defaults('unwarp.estimate.jm');
    realignunwarp.uweoptions.fot        = spm_get_defaults('unwarp.estimate.foe');
    realignunwarp.uweoptions.sot        = spm_get_defaults('unwarp.estimate.soe');
    realignunwarp.uweoptions.uwfwhm     = spm_get_defaults('unwarp.estimate.fwhm');
    realignunwarp.uweoptions.rem        = spm_get_defaults('unwarp.estimate.rem');
    realignunwarp.uweoptions.noi        = spm_get_defaults('unwarp.estimate.noi');
    realignunwarp.uweoptions.expround   = spm_get_defaults('unwarp.estimate.expround');
    realignunwarp.uwroptions.uwwhich    = spm_get_defaults('realign.write.which');
    realignunwarp.uwroptions.rinterp    = spm_get_defaults('realign.write.interp');
    realignunwarp.uwroptions.wrap       = [0 1 0];                                      % wraping along the Y axis (because acquisition in PA direction)
    realignunwarp.uwroptions.mask       = spm_get_defaults('realign.write.mask');
    realignunwarp.uwroptions.prefix     = 'u';
    
    matlabbatch{stage}.spm.spatial.realignunwarp = realignunwarp;
    
end

% OPTION 2: Simple re-alignment
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ismember('realign', actions)
    
    stage = stage + 1;
    clear estwrite
    % get prefix of files, depending on previous processing steps
    if ismember('slicetiming', actions)
        prefix = 'a';
    else
        prefix = '';
    end
    
    for iRun = 1:nrun
        % add prefix to file, depending on previous preprocessing steps
        afflist = AddPrefix(funcfiles, prefix);
        
        % list files for each session
        estwrite.data{iRun} = afflist{iRun};
    end
    
    estwrite.eoptions.quality   = spm_get_defaults('realign.estimate.quality');
    estwrite.eoptions.sep       = spm_get_defaults('realign.estimate.quality');
    estwrite.eoptions.fwhm      = spm_get_defaults('realign.estimate.fwhm');
    estwrite.eoptions.rtm       = spm_get_defaults('realign.estimate.rtm');     % register to 1st
    estwrite.eoptions.interp    = spm_get_defaults('realign.estimate.interp');
    estwrite.eoptions.wrap      = [0 1 0];                                      % wraping along the Y axis (because acquisition in PA direction)
    estwrite.eoptions.weight    = '';
    estwrite.roptions.which     = spm_get_defaults('realign.write.which');
    estwrite.roptions.interp    = spm_get_defaults('realign.write.interp');
    estwrite.roptions.wrap      = [0 1 0];                                      % wraping along the Y axis (because acquisition in PA direction)
    estwrite.roptions.mask      = spm_get_defaults('realign.write.mask');
    estwrite.roptions.prefix    = 'r';
    
    matlabbatch{stage}.spm.spatial.realign.estwrite = estwrite;
    
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
    
    clear preproc
    preproc.channel.vols        = { anatfile };
    preproc.warp.write          = [0 1];
    preproc.channel.biasreg     = 0.001;                    % light regularisation (suggested by default)
    preproc.channel.biasfwhm    = 60;                       % cutoff in mm (suggested by default)
    preproc.channel.write       = [0 1];                    % save bias corrected image
    
    ngaus  = [1 1 2 3 4 2];                                 % (values suggested by default)
    native = [1 1 1 0 0 0];                                 % (values suggested by default)
    for c = 1:6                                             % tissue class c
        preproc.tissue(c).tpm = {
            fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
        preproc.tissue(c).ngaus = ngaus(c);
        preproc.tissue(c).native = [native(c) 0];
        preproc.tissue(c).warped = [0 0];
    end
    
    preproc.warp.mrf            = 1;                        % (suggested by default)
    preproc.warp.cleanup        = 1;                        % (suggested by default)
    preproc.warp.reg            = [0 0.001 0.5 0.05 0.2];   % (suggested by default)
    preproc.warp.affreg         = 'mni';                    % (suggest by default)
    preproc.warp.fwhm           = 0;                        % (suggested by default)
    preproc.warp.samp           = 3;                        % (suggested by default)
    preproc.warp.write          = [1 1];                    % write the forward & inverse deformation fields
    
    matlabbatch{stage}.spm.spatial.preproc = preproc;
    
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
    matlabbatch{stage_Imcal}.spm.util.imcalc.expression     = '(i1 + i2 + i3) .* i4'; % GM, WM, CSF) in skull
    matlabbatch{stage_Imcal}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{stage_Imcal}.spm.util.imcalc.options.dmtx   = 0;    % don't read data as matrix
    matlabbatch{stage_Imcal}.spm.util.imcalc.options.mask   = 0;    % no mask
    matlabbatch{stage_Imcal}.spm.util.imcalc.options.interp = 1;    % tri-linear interpolation    
    matlabbatch{stage_Imcal}.spm.util.imcalc.options.dtype  = 4;    % signed short
    
    
    % coregistration: mean EPI to Brain
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % get the calculated anatomical brain image for this subject
    ImCalBrain = [adir, 'Brain.nii,1'];
    
    % get EPIs
    % List the EPI, with ua prefix
    prefix = '';
    if  ismember('slicetiming', actions) && ismember('realign', actions) &&  ismember('topup', actions); prefix = 'tra'; end
    if ~ismember('slicetiming', actions) && ismember('realign', actions) &&  ismember('topup', actions); prefix = 'tr';  end
    if  ismember('slicetiming', actions) && ismember('unwrap', actions)  &&  ismember('topup', actions); prefix = 'tua'; end
    if ~ismember('slicetiming', actions) && ismember('unwrap', actions)  &&  ismember('topup', actions); prefix = 'tu';  end
    if  ismember('slicetiming', actions) && ismember('realign', actions) && ~ismember('topup', actions); prefix = 'ra';  end
    if ~ismember('slicetiming', actions) && ismember('realign', actions) && ~ismember('topup', actions); prefix = 'r';   end
    if  ismember('slicetiming', actions) && ismember('unwrap', actions)  && ~ismember('topup', actions); prefix = 'ua';  end
    if ~ismember('slicetiming', actions) && ismember('unwrap', actions)  && ~ismember('topup', actions); prefix = 'u';   end
    OtherEPI = AddPrefix(funcfiles, prefix);
    
    % compute mean image
    stage = stage + 1;
    matlabbatch{stage}.spm.util.imcalc.input            = vertcat(OtherEPI{:});
    matlabbatch{stage}.spm.util.imcalc.output           = 'mean_preproc_EPI';
    matlabbatch{stage}.spm.util.imcalc.outdir           = {fdir};
    matlabbatch{stage}.spm.util.imcalc.expression       = 'mean(X)';
    matlabbatch{stage}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
    matlabbatch{stage}.spm.util.imcalc.options.dmtx     = 1; % read images as matrix
    matlabbatch{stage}.spm.util.imcalc.options.mask     = 0; % no masking
    matlabbatch{stage}.spm.util.imcalc.options.interp   = 1;
    matlabbatch{stage}.spm.util.imcalc.options.dtype    = 4;
    
    stage = stage + 1;
    OtherEPI = AddPrefix(cffiles, prefix);
    clear coreg
    coreg.estimate.ref(1)            = {ImCalBrain};
    coreg.estimate.source(1)         = {[fdir, 'mean_preproc_EPI.nii,1']};
    coreg.estimate.other             = OtherEPI;
    coreg.estimate.eoptions.cost_fun = spm_get_defaults('coreg.estimate.cost_fun');
    coreg.estimate.eoptions.sep      = spm_get_defaults('coreg.estimate.sep');
    coreg.estimate.eoptions.tol      = spm_get_defaults('coreg.estimate.tol');
    coreg.estimate.eoptions.fwhm     = spm_get_defaults('coreg.estimate.fwhm');
    
    matlabbatch{stage}.spm.spatial.coreg = coreg;
    
    
    % Spatial normalisation of EPIs, given the normalized segmentation (and its
    % affine transformation)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % get (previously computed) deformation field
    ForwardDef = [adir, 'y_', spm_select('List', adir, regexp_anat)];
    
    stage = stage + 1;
    
    matlabbatch{stage}.spm.spatial.normalise.write.subj.def         = {ForwardDef};
    matlabbatch{stage}.spm.spatial.normalise.write.subj.resample    = OtherEPI;
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb      = spm_get_defaults('normalise.write.bb');
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.vox     = voxel_size;
    matlabbatch{stage}.spm.spatial.normalise.write.woptions.interp  = spm_get_defaults('normalise.write.interp');
    
    
    % Spatial normalization of anat, given the normalized segmentation (and its
    % affine transformation)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
    
    matlabbatch{stage}.spm.spatial.smooth.data   = NormEPIs;
    matlabbatch{stage}.spm.spatial.smooth.fwhm   = smoothing_kernel;
    matlabbatch{stage}.spm.spatial.smooth.dtype  = 0; % output data type = same as input
    matlabbatch{stage}.spm.spatial.smooth.im     = 0; % no implicit masking
    matlabbatch{stage}.spm.spatial.smooth.prefix = 's';
end

% =========================================================================
%                            SAVE SPM BATCH
% =========================================================================
matfile = sprintf('%s/batch_preproc_sub%02.0f.mat', datadir, iSub);
if exist(matfile, 'file')
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
