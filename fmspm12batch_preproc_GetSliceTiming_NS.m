function [SliceTiming, TR, TE, SliceThickness, SpacingBetweenSlices, ...
            NumberOfSlices, PixelSpacing, total_readout_time_spm, total_readout_time_fsl] ...
            = fmspm12batch_preproc_GetSliceTiming_NS(iSub, funcdir, datadir, spm_path, ForceToThisDICOM)
% Function to get (and save) the slice timing info from the Neurospin
% server for a particular subject. The timing is read in the header of a
% parciular dicom volume; it should not vary significantly across volume
% and subject (the accuracy of these numbers is <5ms).

% The 5th argument is optinal. If provided, it should be the full path name
% of the DICOM file serving as reference.

% add spm in the path
addpath(spm_path)

% get the directory of imported data for this subject
subjdir = sprintf('%s/subj%02.0f/%s/', datadir, iSub, funcdir);

if nargin < 5
    
    % Find NIP and date of aquisition
    % Assume that all file names follow the pattern *_NIP_subDate_*.nii,
    % with NIP made of 2 letters and 6 digits, and subDate made of 8 digits.
    % File names have this form when they are imported with Christophe
    % Pallier's script.
    fname = spm_select('List', subjdir, '^*\.nii');
    fname = deblank(fname(1,:));
    ind = regexp(fname, '\_[a-z]{2}\d{6}\_');
    if isempty(ind) || numel(ind) > 1; error('cannot find NIP'); end;
    NIP = fname(ind+1:ind+8);
    subDate = fname(ind+10:ind+17);
    
    % go into folder on the neurospin server
    if str2double(subDate(1:4)) > 2015 || (str2double(subDate(1:4)) == 2015 && str2double(subDate(1:4))>= 11)
        % data recorder after 2015/11 are on the Prisma repository
        base_ns = ['/neurospin/acquisition/database/Prisma_fit/', subDate];
    else
        % data recorder before 2015/11 are on the Trio repository
        base_ns = ['/neurospin/acquisition/database/TrioTim/', subDate];
    end
    dname = dir([base_ns, '/', NIP, '*']);
    dname = dname.name;
    subjdir_ns = [base_ns, '/', dname];
    
    % get name of the the 1st mbepi folder (that is not SBref!!)
    dname = dir(subjdir_ns);
    flag = 0;
    for iDir = 1:length(dname)
        if ~isempty(strfind(dname(iDir).name, 'mbepi')) && isempty(strfind(dname(iDir).name, 'SBRef'))
            flag = 1;
            break
        end
    end
    if flag == 0
        error('No directory found for mbepi!!')
    end
    sessname = dname(iDir).name;
    
    % Load the 2nd dicom header (not the 1st, as recommended by the SPM manual
    % on slice timing correction)
    sess_dir_ns = [subjdir_ns, '/', sessname, '/'];
    flist = dir([sess_dir_ns, '*.dcm']);
    
    hdr = spm_dicom_headers([sess_dir_ns, flist(2).name]);
    
else
    % FORCE TO TAKE THIS DICOM INTO ACCOUNT
    hdr = spm_dicom_headers(ForceToThisDICOM);
end

% Get the timing info:
SliceTiming = hdr{1}.Private_0019_1029; % slice timing
TR = hdr{1}.RepetitionTime/1000; % Convert to seconds
TE = hdr{1}.EchoTime;            % leave in ms.
SliceThickness = hdr{1}.SliceThickness;
SpacingBetweenSlices = hdr{1}.SpacingBetweenSlices;
NumberOfSlices = hdr{1}.Private_0019_100a;
PixelSpacing = hdr{1}.PixelSpacing; % [x, y] resolution in mm

% GET TOTAL READOUT TIME
% ----------------------
% explanation can be found in several web pages:
% https://lcni.uoregon.edu/kb-articles/kb-0003
% http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TOPUP/TopupUsersGuide
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;92bd6f89.1403

% BW is the effective (or "reconstructed") BW in PE direction
% nVx is the reconstructed number of line
% 1/(BW*nVx) is the effective echo time. It is indep. from the GRAPPA 
% since both BW and nVx are in the reconstructed space.
% The echo time reported by Siemens is the actual echo time. The actual
% time corresponds to the effective time times the GRAPPA factor. 
% To check that the computation is correct, the effective echo time can be 
% compared with the actual echo time interval reported in the Siemens PDF. 
% actual echo time (Siemens) = iPAT/(BW*nVx).
%
% FSL and SPM needs is the effective echo time with a minor difference:
%   SPM counts from start of first echo to end of last echo (= exactly the echo time)
%   FLS counts from the middel of first echo to end of last echo.
%
% NB: FSL (topup guide) say that we should count with the number of
% "reconstructed" echos, not the actual one (in other words, which should
% not care about the GRAPPA factor).

% get number of voxel in the encoding direction
nVx = hdr{1}.Private_0051_100b;
nVx = str2double(nVx(1:(strfind(nVx, '*')-1)));

% get "Bandwidth per pixel phase encode" in Hz
BW = hdr{1}.Private_0019_1028;

total_readout_time_spm = 1/BW;
total_readout_time_fsl = (nVx-1)/(BW*nVx);

% save the sclice timing info
fname = [subjdir, 'SliceTimingInfo.mat'];
if exist(fname)
    delete(fname); fprintf('\nremove previous SliceTimingInfo.mat\n')
end
save(fname, 'SliceTiming', 'TR', 'TE', 'SliceThickness', 'SpacingBetweenSlices', ...
    'NumberOfSlices', 'PixelSpacing', 'total_readout_time_spm', 'total_readout_time_fsl')

