function [SliceTiming, TR, TE, SliceThickness, SpacingBetweenSlices, NumberOfSlices] ...
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

% get subject NIP and date
fname = spm_select('List', subjdir, '^B0_1_.*\.nii');
ind_ = strfind(fname, '_');
subID = fname(ind_(2)+1:ind_(3)-1);
subDate = fname(ind_(3)+1:ind_(4)-1);

% go into folder on the neurospin server
if str2double(subDate(1:4)) > 2015 || (str2double(subDate(1:4)) == 2015 && str2double(subDate(1:4))>= 11)
    % data recorder after 2015/11 are on the Prisma repository
    base_ns = ['/neurospin/acquisition/database/Prisma_fit/', subDate];
else
    % data recorder before 2015/11 are on the Trio repository
    base_ns = ['/neurospin/acquisition/database/TrioTim/', subDate];
end
dname = dir([base_ns, '/', subID, '*']);
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

% save the sclice timing info
fname = [subjdir, 'SliceTimingInfo.mat'];
if exist(fname) 
    delete(fname); fprintf('\nremove previous SliceTimingInfo.mat\n')
end
save(fname, 'SliceTiming', 'TR', 'TE', 'SliceThickness', 'SpacingBetweenSlices', 'NumberOfSlices')

