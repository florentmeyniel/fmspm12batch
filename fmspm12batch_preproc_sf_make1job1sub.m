function [matfile, matlabbatch] = fmspm12batch_preproc_sf_make1job1sub(iSub, dataloc, param)

% TO IMPROVE: l. 166 : faire une boucle sur les session pour récupérer la première image automatiquement!!!

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

% unpack input data
% =========================================================================
spm_path         = dataloc.spm_path;            % spm directory  
datadir          = dataloc.datadir;             % root directory for subject data
regexp_func      = dataloc.regexp_func;         % regular expression to recognize functional sessions to analyze 
regexp_anat      = dataloc.regexp_anat;         % regular expression to recognize T1
funcdir          = dataloc.funcdir;             % path of fMRI data (4D nifti) within subject directory
anatdir          = dataloc.anatdir;             % path of anatomical image within subject directory

nslices          = param.nslices;               % number of slices
deltaEPI         = param.deltaEPI;              % readout between 2 EPI in ms ('Ecart echo' in the Siemens PDF)
iPAT             = param.iPAT;                  % EPI acceleration
voxel_size       = param.voxel_size; 
smoothing_kernel = param.smoothing_kernel;      % 1st level smoothing
TR               = param.TR;
slice_timing     = param.slice_timing;

% initialize spm
% =========================================================================
addpath(spm_path)
spm('defaults', 'FMRI');
spm_jobman('initcfg');
stage = 0;

% Get the data files
% =========================================================================

% Get subject directory (with anat & fMRI)
subjdir = sprintf('%s/subj%02.0f/', datadir, iSub);

% Get T1 file to analyze
adir =  sprintf('%s/%s/', subjdir, anatdir);
anatfile = spm_select('FPList', adir, regexp_anat);
if isequal(anatfile,  '')
    warning(sprintf('No anat file found for %s', ...
        subjects{csubj}))
    return
end

% Get functional files to analyze
fdir =  sprintf('%s%s/', subjdir, funcdir);
ffiles = spm_select('List', fdir, regexp_func);
nrun = size(ffiles,1);
if nrun == 0
    warning(sprintf('No functional file found for %s', ...
        subjects{csubj}))
    return
end
funcfiles = cell(nrun, 1);
cffiles = cellstr(ffiles);
for i = 1:nrun
    funcfiles{i} = cellstr(spm_select('ExtFPList', fdir, ['^', cffiles{i}], Inf));
end


% Slice timing
% =========================================================================
stage = stage + 1;
stage_slicetiming = stage;

matlabbatch{1}.spm.temporal.st.scans = { funcfiles{:} };
matlabbatch{1}.spm.temporal.st.nslices = nslices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = 0;              % TA can be left to 0 if slice timing is explit (as ms, not order)
matlabbatch{1}.spm.temporal.st.so = slice_timing;
matlabbatch{1}.spm.temporal.st.refslice = 0;        % reference slice is actually ref. time (in ms) because slice times are provided.
matlabbatch{1}.spm.temporal.st.prefix = 'a';


% Compute FieldMap correction
% =========================================================================
stage = stage + 1;
stage_fieldmap = stage;

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
cmd = sprintf('cd %s; fsl5.0-fslroi %s %s 0 1', ...
    [fdir, 'B0'], ...
    fname_B0mag, strcat([fname_B0mag(1:7), 'shortTE_', fname_B0mag(8:end)])); 
unix(cmd)
fname_B0phase = spm_select('List', [fdir, 'B0'], '^B0_phase.*\.nii');
cmd = sprintf('cd %s; fsl5.0-fslroi %s %s 0 1', ...
    [fdir, 'B0'], ...
    fname_B0phase, strcat([fname_B0phase(1:9), 'shortTE_', fname_B0phase(10:end)])); 
unix(cmd)

% 2nd (1) = intermediate: 5.46 ms
cmd = sprintf('cd %s; fsl5.0-fslroi %s %s 1 1', ...
     [fdir, 'B0'], ...
    fname_B0mag, strcat([fname_B0mag(1:7), 'longTE_', fname_B0mag(8:end)])); 
unix(cmd)
cmd = sprintf('cd %s; fsl5.0-fslroi %s %s 1 1', ...
     [fdir, 'B0'], ...
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
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.et = [3 5.46];                     % [shortTE longTE] in ms
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.blipdir = 1;                       % +1 for P -> A (would be -1 for A -> P)
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.tert = nslices*deltaEPI/iPAT;      % EPI readout time
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.epifm = 0;                         % 0: fieldmap is not an EPI image (it is a gradiant echo)
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.ajm = 0;                           % no jacobian modulation (as recommanded by SPM)
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.template = ...
                                                                {[spm_path, '/toolbox/FieldMap/T1.nii']};   % for visual comparison 
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.session(1).epi = {funcfiles{1}{1}};                     % specify 1st EPI for session alignement (for quality control with a display)
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.session(2).epi = {funcfiles{2}{1}};                     % specify 1st EPI for session alignement (for quality control with a display)
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.session(3).epi = {funcfiles{3}{1}};                     % specify 1st EPI for session alignement (for quality control with a display)
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.session(4).epi = {funcfiles{4}{1}};                     % specify 1st EPI for session alignement (for quality control with a display)
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.matchvdm = 1;                                           % so that the VDM is aligned on the 1st EPI of each session
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.sessname = 'session';                                   % VDM file extension
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.writeunwarped = 1;
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.anat = {[anatfile, ',1']};
matlabbatch{stage}.spm.tools.fieldmap.phasemag.subj.matchanat = 1;                                          % for visual comparison


% Realign and unwrap the functional images
% =========================================================================
% Note that in SPM, if multiple sessions are provided, the 1st scan of each 
% session are realigned to each other (on the first scan of the first session)
% and then, within each session, the scan are realigned to the 1st scan of that 
% session. The realignment parameters saved are thus relative to the 1st scan
% of the 1st session, which is not a problem since the 1st scan of the different
% sessions are first aligned to each other.
stage = stage + 1;
stage_realign = stage;

for iRun = 1:nrun
    matlabbatch{stage}.spm.spatial.realignunwarp.data(iRun).scans(1) = ...
        cfg_dep(sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)', iRun), ...
        substruct('.','val', '{}',{stage_slicetiming}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('()',{iRun}, '.','files'));
    matlabbatch{stage}.spm.spatial.realignunwarp.data(iRun).pmscan(1) = ...
        cfg_dep(sprintf('Phase and Magnitude Data: Voxel displacement map (Subj 1, Session %d)', iRun), ...
        substruct('.','val', '{}',{stage_fieldmap}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('()',{1}, '.','vdmfile', '{}',{iRun}));    
end

matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.ewrap = [0 1 0];                                      % wraping along the Y axis (because acquisition in PA direction)
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.wrap = [0 1 0];
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

% store dependencies in a structure: the realigned & unwrapped data will be
% used for the co-registration step and the normalization step.
clear sessdep
for iRun = 1:nrun
    sessdep(iRun) = ...
        cfg_dep(sprintf('Realign & Unwarp: Unwarped Images (Sess %d)', iRun), ...
        substruct('.','val', '{}',{stage_realign}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','sess', '()',{iRun}, '.','uwrfiles'));
end

% NB: about TOPUP (FSL) corrections
% Say that function images are acquired in the P -> A orientation, and that
% one has two calibration files PA and AP.
% 1) the PA file should be aligned with the functional file. This can be
% achieved by adding the PA as an extra session (containing only one image)
% to the previous realign & unwrap step (note: it is the realign step that
% is needed, unwrap is for another purpose - the FieldMap correction).
% 2) compute deformation map with topup on AP & PA
% 3) apply the deformation map to sessions (applytopup)


% Segmentation and spatially normalize anatomy with an MNI template by
% computing an affine transformation
% =========================================================================
% (using SPM default parameters)

stage = stage + 1;
stage_segmentation = stage;

matlabbatch{stage}.spm.spatial.preproc.channel.vols = { anatfile };
matlabbatch{stage}.spm.spatial.preproc.warp.write = [0 1];
matlabbatch{stage}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{stage}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{stage}.spm.spatial.preproc.channel.write = [0 1];
ngaus  = [1 1 2 3 4 2];
native = [1 1 1 0 0 0];
for c = 1:6 % tissue class c
    matlabbatch{stage}.spm.spatial.preproc.tissue(c).tpm = {
        fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
    matlabbatch{stage}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
    matlabbatch{stage}.spm.spatial.preproc.tissue(c).native = [native(c) 0];
    matlabbatch{stage}.spm.spatial.preproc.tissue(c).warped = [0 0];
end
matlabbatch{stage}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{stage}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{stage}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{stage}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{stage}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{stage}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{stage}.spm.spatial.preproc.warp.write = [1 1];

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


% coregistration: mean function to Brain
% =========================================================================
% (using SPM default parameters)

stage = stage + 1;
stage_coregister = stage;

matlabbatch{stage_coregister}.spm.spatial.coreg.estimate.ref(1) = ...
    cfg_dep('Image Calculator: Imcalc Computed Image', ...
    substruct('.','val', '{}',{stage_Imcal}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
matlabbatch{stage_coregister}.spm.spatial.coreg.estimate.source(1) = ...
    cfg_dep('Realign & Unwarp: Unwarped Mean Image', ...
    substruct('.','val', '{}',{stage_realign}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','meanuwr'));

matlabbatch{stage_coregister}.spm.spatial.coreg.estimate.other = sessdep;
matlabbatch{stage_coregister}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{stage_coregister}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{stage_coregister}.spm.spatial.coreg.estimate.eoptions.tol = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{stage_coregister}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];


% Spatial normalisation of EPIs, given the normalized segmentation (and its
% affine transformation)
% =========================================================================
% (using SPM default parameters)

stage = stage + 1;
stage_normalisationEPI = stage;

matlabbatch{stage}.spm.spatial.normalise.write.subj.def(1) = ...
    cfg_dep('Segment: Forward Deformations', ...
    substruct('.','val', '{}',{stage_segmentation}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','fordef', '()',{':'}));
matlabbatch{stage}.spm.spatial.normalise.write.subj.resample = sessdep;
matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
matlabbatch{stage}.spm.spatial.normalise.write.woptions.vox = voxel_size;
matlabbatch{stage}.spm.spatial.normalise.write.woptions.interp = 4;


% Spatial normalization of anat, given the normalized segmentation (and its
% affine transformation)
% =========================================================================
% (using SPM default parameters)

stage = stage + 1;
stage_normalisationAnat = stage;

matlabbatch{stage}.spm.spatial.normalise.write.subj.def(1) = ... 
    cfg_dep('Segment: Forward Deformations', ...
    substruct('.','val', '{}',{stage_segmentation}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','fordef', '()',{':'}));
matlabbatch{stage}.spm.spatial.normalise.write.subj.resample(1) = ...
    cfg_dep('Image Calculator: Imcalc Computed Image', ...
    substruct('.','val', '{}',{stage_Imcal}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','files'));
matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
matlabbatch{stage}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{stage}.spm.spatial.normalise.write.woptions.interp = 4;


% Smoothing
% =========================================================================
stage = stage + 1;
stage_smoothing = stage;

matlabbatch{stage}.spm.spatial.smooth.data(1) = ...
    cfg_dep('Normalise: Write: Normalised Images (Subj 1)', ...
    substruct('.','val', '{}',{stage_normalisationEPI}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('()',{1}, '.','files'));
matlabbatch{stage}.spm.spatial.smooth.fwhm = smoothing_kernel;
matlabbatch{stage}.spm.spatial.smooth.dtype = 0;
matlabbatch{stage}.spm.spatial.smooth.im = 0;
matlabbatch{stage}.spm.spatial.smooth.prefix = 's';


% Save SPM batch
% =========================================================================
matfile = sprintf('%s/batch_preproc_sub%02.0f.mat', datadir, iSub);
if exist(matfile) 
    delete(matfile); fprintf('\nremove previous batch\n')
end
fprintf('\n Writing SPM batch for subject %d:', iSub)
fprintf('\n %s\n', matfile)
save(matfile,'matlabbatch');

