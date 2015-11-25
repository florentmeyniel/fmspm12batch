function fmspm12batch_CheckROI(ClusterCoord_vx, ROIname, vx2mni)
% Function to plot the ROIs
% fmspm12batch_CheckROI(ClusterCoord_vx, ROIname, vx2mni)
%   ClusterCoord_vx: cell of 3 x n voxel coordinate (voxel space)
%   ROIname: cell of ROI names
%   vx2mni: spm MNI 2 voxel 4x4 matrix


if iscell(ClusterCoord_vx)
    nROI = numel(ClusterCoord_vx);
    for iROI = 1:nROI
        % convert to mni
        XYZ = ClusterCoord_vx{iROI}';
        mni = vx2mni*[XYZ(:,1) XYZ(:,2) XYZ(:,3) ones(size(XYZ,1),1)]';
        mni = mni';
        mni(:,4) = [];
        XYZ = mni;
        
        figure; set(gcf, 'Color', [1 1 1])
        set(gcf, 'Name', ROIname{iROI})
        spm_mip(ones(size(XYZ,1), 1), XYZ, vx2mni);
        colormap('bone')
    end
else
    % convert to mni
    XYZ = ClusterCoord_vx';
    mni = vx2mni*[XYZ(:,1) XYZ(:,2) XYZ(:,3) ones(size(XYZ,1),1)]';
    mni = mni';
    mni(:,4) = [];
    XYZ = mni;
    
    figure; set(gcf, 'Color', [1 1 1])
    set(gcf, 'Name', ROIname)
    spm_mip(ones(size(XYZ,1), 1), XYZ, vx2mni);
    colormap('bone')
end
