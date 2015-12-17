function [ClusterCoord_vx, xSPM] = fmspm12batch_GetClusterCoord(Peakmm, xSPM)
% Function to compute the list of voxels (in the voxel space) of the
% cluster corresponding to the input coordinate (in mm), and in the SPM.mat
% selected interactively with its statistical threshold.
% NB: Take the coordinate from the table, or from the Result panel, but not
% the local maximal printed in the command line (which does not correspond
% to the actual grid of voxels).
%
% Usage:
% [ClusterCoord_vx, xSPM] = fmspm12batch_GetClusterCoord(Peakmm)
% NB: an extra argument xSPM can be passed in:
%   - it can the output from a previous call of fmspm12batch_GetClusterCoord
%     to avoid entering again which model, stat, ... should be used
%     (typically, if multiple ROI are defined for the same model)
%   - it can be a xSPM with the following field (easy to batch):
%       .swd: full path of the SPM.mat
%       .Ic: number of the contrast (cf. ordering of spmT maps)
%       .k: minimum extend of clusters
%       .u p-value to threshold the map
%       .thresDesc: type of stat: 'FWE', 'FDR' or 'none'

% initialize
try ver = spm('Version');
catch
    error('the SPM toolbox is not in the Matlab''s path')
end
spm_results_ui('clear')
% spm_results_ui('close')


% Get the xSPM, let the user select the contrast and statistical threshold.
if nargin == 1
    [~, xSPM] = spm_getSPM;
else
    if ~isfield(xSPM, 'title')
        % this means that xSPM was not computed automatically by SPM,
        % but instead provided with minimal information by the user.
        
        % complete some fields with default settings
        xSPM.n  = 1;  % value >1 are used to compute conjunctions
        xSPM.Im = []; % empty because there is no masking
        xSPM.pm = []; % no masking = no p-value to define the mask
        xSPM.Ex = []; % irrelevant flag when there is no masking
        
        % get the contrast title field
        SPM = load([xSPM.swd, '/SPM.mat']);
        SPM = SPM.SPM;
        xSPM.title = SPM.xCon(xSPM.Ic).name;
        
        % let SPM finish the job!
        [~, xSPM] = spm_getSPM(xSPM);
    end
end

% List all the cluster in this thresholded map
clusters = spm_clusters(xSPM.XYZ);

% Get the voxel, in the list of coordinate, corresponding to a peak.
% NB: Take the coordinate from the table, or from the Result panel, but not
% the local maximal printed in the command line (which does not correspond
% to the actual grid of voxels).
peakloc = ismember(xSPM.XYZmm', Peakmm,'rows');

% check that there is a unique location
if sum(peakloc) ~= 1
    error('the location correspond to %d (not 1) voxel', sum(peakloc))
end

% Get coordinate of all voxels in this cluster
ClusterCoord_vx = xSPM.XYZ(:, clusters==clusters(peakloc));