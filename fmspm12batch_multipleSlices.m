% This script plots the same anatomical slice, and overlay the results from
% different subjects.
% It prompts the parameters to set (files, which orientation, which slice,
% ...).
% It sets up a slover structure (as used in SPM) and call the modified
% slover function paint_fm.m.


% initialize
try ver = spm('Version'); 
catch
    error('the SPM toolbox is not in the Matlab''s path')
end
spm_results_ui('clear')
spm_results_ui('close')

%-Get individual SPMs
%--------------------------------------------------------------------------

% get the rfx model directory (to get the model name)
modeldir = spm_select(1,'dir','Select rfx model');
ind = strfind(modeldir, '/');
modelname = modeldir(ind(end)+1:end);
ind = strfind(modeldir, '/group_analysis/');
rootdir = modeldir(1:ind-1);

% get list of subjects to plot
sublist = spm_input('subject list','+1','e','');
nSub = numel(sublist);

% built list of individual SPM.mat
SPMs = [];
for iSub = 1:nSub
    SPMs(iSub,:) = ...
        sprintf('%s/subj%02.0f/first_level_estimates/%s/SPM.mat', ...
        rootdir, sublist(iSub), modelname);
end
SPMs = char(SPMs);


%-choose anat file for rendering
%--------------------------------------------------------------------------
% initialize the slover structure
so = slover;

anat_auto = spm_input('Image for rendering on', '+1', ...
    'Auto|Manual', [1, 0], 1);

if anat_auto
    img = sprintf('%s/group_analysis/template/MeanT1.nii', rootdir);
    if ~exist(img, 'file')
        [img,sts] = spm_select(1,'image','Can'' find anat file automatically. Please select');
    end
else
    [img,sts] = spm_select(1,'image','Select image for rendering on');
end

% enter the rendering file in the slover structure
so.img.vol = spm_vol(img);
so.img.prop = 1;


%-Get transformation (which plane to plot)
%--------------------------------------------------------------------------
so.transform = deblank(spm_input('Image orientation', '+1', ...
    'Axial|Coronal|Sagittal', char('axial','coronal','sagittal'), 1));
so = fill_defaults(so);
slices = so.slices;
so.slices = spm_input('Slices to display (mm)', '+1', 'e', ...
    sprintf('%0.0f:%0.0f:%0.0f',slices(1),mean(diff(slices)),slices(end)));

so.figure = spm_figure('GetWin', 'SliceOverlay');


%-Get input parameter xSPM (1st level results)
%--------------------------------------------------------------------------
SPMs = cellstr(SPMs);
clear xSPM
for iSub = 1:nSub
    if iSub == 1
        % take 1st subject as reference. This prompts for a statistical threshold.
        xSPM.swd = spm_file(SPMs{iSub},'fpath');
        [tmp, xSPM] = spm_getSPM(xSPM);
        
        % get the statistical threshold (like in spm_check_results.m)
        td = regexp(xSPM.thresDesc,'p\D?(?<u>[\.\d]+) \((?<thresDesc>\S+)\)','names');
        if isempty(td)
            td = regexp(thd,'\w=(?<u>[\.\d]+)','names');
            td.thresDesc = 'none';
        end
        if strcmp(td.thresDesc,'unc.'), td.thresDesc = 'none'; end
        xSPM.u          = str2double(td.u);
        xSPM.thresDesc  = td.thresDesc;

        % add this participant to the 'slover' object
        so = add_spm(so, xSPM);
    else
        % estimation of the subject SPM.mat with the statistical threshold
        % defined for the 1st subject
        xSPM.swd = spm_file(SPMs{iSub},'fpath');
        [tmp, myxSPM] = spm_getSPM(xSPM);
        
        % add this participant to the 'slover' object
        so = add_spm(so, myxSPM);
    end
end


%-Assign the same color scale to all subjects
%--------------------------------------------------------------------------
% get min & max over subjects
mini = so.img(2).range(1); maxi = so.img(2).range(2);
for iSub = 2:nSub;
    mini = min([mini, so.img(1+iSub).range(1)]);
    maxi = max([maxi, so.img(1+iSub).range(2)]);
end

% set this max for all subjects
for iSub = 1:nSub;
    so.img(1+iSub).range(1) = mini;
    so.img(1+iSub).range(2) = maxi;
end

% trick the paint_fm function, predenting that there is only 1 colorar bar,
% that of figure 2, to plot activations
so.cbar = 2;


%-Plot the results
%--------------------------------------------------------------------------
so = paint_fm(so);

