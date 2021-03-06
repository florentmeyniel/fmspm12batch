function fmspm12batch_preproc_newArchi(paramfile)
% Function to copy the smoothed, normalized files into a new directory.
% Load the preprocessing parameter file and assume that the
% subject-specific folder are immediately in a folder labelled 'raw_data'.
% 
% syntax: fmspm12batch_preproc_newArchi(paramfile)

% Get the parameters
eval(paramfile);

% Make the directory where data to analyzed are located
ind = strfind(datadir, 'raw_data');
newdirname = [datadir(1:ind-1), 'analyzed_data'];
if ~exist(newdirname); mkdir(newdirname); end

for iSub = 1:length(sublist)
    
    % funtional files
    % ===============
    
    % list preprocessed files to move (sw*.nii)
    orign_subjdir = sprintf('%s/subj%02.0f/%s/', datadir, sublist(iSub), 'fMRI');
    flist = cellstr(spm_select('List', orign_subjdir, '^swt.*\.nii'));
    
    % make a new directory
    target_subjdir = sprintf('%s/subj%02.0f/%s/', newdirname, sublist(iSub), 'preprocEPI');
    if ~exist(target_subjdir, 'dir')
        mkdir(target_subjdir)
    end
    
    % move files
    for iFile = 1:length(flist)
        movefile([orign_subjdir, flist{iFile}], [target_subjdir, flist{iFile}])
    end
    
    % side information: slice timing & movement parameters
    % ====================================================
    % Slice timing info
    copyfile([orign_subjdir, 'SliceTimingInfo.mat'], [target_subjdir, 'SliceTimingInfo.mat'])
    
    movparfile = cellstr(spm_select('List', orign_subjdir, '^rp_aepi_.*\.txt'));
    for iFile = 1:length(movparfile)
        copyfile([orign_subjdir, movparfile{iFile}], ...
            [target_subjdir, movparfile{iFile}])
    end
        
    % anat files
    % ==========
    % get file to move (sw*.nii)
    orign_subjdir = sprintf('%s/subj%02.0f/%s/', datadir, sublist(iSub), 'anat');
    fname = 'wBrain.nii';
    
    % make a new directory
    target_subjdir = sprintf('%s/subj%02.0f/%s/', newdirname, sublist(iSub), 'preprocAnat');
    if ~exist(target_subjdir, 'dir')
        mkdir(target_subjdir)
    end
    
    % copy files
    copyfile([orign_subjdir, fname], [target_subjdir, fname])
    
end
