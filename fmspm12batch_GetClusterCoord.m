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