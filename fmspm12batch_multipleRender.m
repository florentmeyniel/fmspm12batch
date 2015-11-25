function varargout = fmspm12batch_multipleRender
% Render blobs on a cortical surface.
% This function is modified from spm_render.m to plot a series of 
% individual results
%
% MODIFICATIONS BY FM:
%   => the renderer file is loaded automatically. If not, it is prompted.
%   => the renderer view is prompted
%   => the function automatically search the folder to find the
%   subjects (the list of subject is prompted). The folder architecture
%   must follow that the fmspm12batch toolbox.
%   => the individuals SPM are overlaid on the same rendering view
%

% initialize the SPM window.
try ver = spm('Version'); 
catch
    error('the SPM toolbox is not in the Matlab''s path')
end
spm_results_ui('clear')
spm_results_ui('close')


SVNrev = '$Rev: 6190 $';

global prevrend
if ~isstruct(prevrend)
    prevrend = struct('rendfile','', 'brt',[], 'col',[]);
end

varargout = {};

%-Parse arguments, get data if not passed as parameters
%==========================================================================

%-Get data
%--------------------------------------------------------------------------

% Select the rendering view
rendtype = spm_input('Image orientation', '+1', ...
    'T|Bot|L|R|F|Back', 1:6, 1);

%-Get individual SPMs
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

%- Get stat results for each subject
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
        
        % Add this participant to the data set
        dat(iSub) = struct('XYZ', xSPM.XYZ,...
            't',   xSPM.Z',...
            'mat', xSPM.M,...
            'dim', xSPM.DIM);
    else
        % estimation of the subject SPM.mat with the statistical threshold
        % defined for the 1st subject
        xSPM.swd = spm_file(SPMs{iSub},'fpath');
        [tmp, myxSPM] = spm_getSPM(xSPM);
        
        % Add this participant to the data set
        dat(iSub) = struct('XYZ', myxSPM.XYZ,...
            't',   myxSPM.Z',...
            'mat', myxSPM.M,...
            'dim', myxSPM.DIM);
    end
end


%-Get surface
%--------------------------------------------------------------------------
% try to get the surface automatically
rendfile = [];
try
    ind = strfind(modeldir, 'rfxmodel');
    if ~isempty(ind)
        putative_renderer = dir(strcat([modeldir(1:ind-1), '/template/render*.mat']));
        if numel(putative_renderer) == 1
            rendfile = sprintf('%stemplate/%s', modeldir(1:ind-1), putative_renderer.name);
        end
    else
        ind = findstr(xSPM.swd, '/');
        putative_renderer = dir(strcat([modeldir(1:ind(end-1)), 'preprocAnat/render*.mat']));
        if numel(putative_renderer) == 1
            rendfile = sprintf('%spreprocAnat/%s', xSPM.swd(1:ind(end-1)), putative_renderer.name);
        end
    end
end
if isempty(rendfile)
    spm('FnBanner',mfilename,SVNrev);
    [rendfile, sts] = spm_select(1,{'mat','mesh'},'Render file');
    if ~sts, return; end
end
prevrend.rendfile = rendfile;

ext = spm_file(rendfile,'ext');
loadgifti = false;
if strcmpi(ext,'mat')
    load(rendfile);
    if ~exist('rend','var') && ~exist('Matrixes','var')
        loadgifti = true;
    end
end
if ~strcmpi(ext,'mat') || loadgifti
    try
        rend = gifti(rendfile);
    catch
        error('Cannot read  render file "%s".\n', rendfile);
    end
    if ~isfield(rend,'vertices')
        try
            M = rend;
            rend = gifti(rend.private.metadata(1).value);
            try, rend.cdata = M.cdata(); end
        catch
            error('Cannot find a surface mesh to be displayed.');
        end
    end
    rend = export(rend,'patch');
end

% select rendering view.
rend = rend(rendtype);

%-Get brightness & colours (image)
%--------------------------------------------------------------------------
if nargin < 2  || isempty(prevrend.brt)
    brt = 1;
    brt = spm_input('Style',1,'new|old',[1 NaN], 1);
    
    if isfinite(brt)
        brt = spm_input('Brighten blobs',1,'none|slightly|more|lots',[1 0.75 0.5 0.25], 1);
        col = eye(3);
        % ask for custom colours & get rgb values
        %------------------------------------------------------------------
        if spm_input('Which colours?','!+1','b',{'RGB','Custom'},[0 1],1)
            for k = 1:num
                col(k,:) = uisetcolor(col(k,:),sprintf('Colour of blob set %d',k));
            end
        end
    else
        col = [];
    end
elseif isfinite(brt) && isempty(prevrend.col)
    col = eye(3);
elseif isfinite(brt)  % don't need to check prevrend.col again
    col = prevrend.col;
else
    col = [];
end
prevrend.brt = brt;
prevrend.col = col;

%-Perform the rendering
%==========================================================================
spm('Pointer','Watch');
showbar = 1;
if showbar, spm_progress_bar('Init', size(dat,1)*length(rend),...
            'Formatting Renderings', 'Number completed'); end
for i=1:length(rend)
    rend{i}.max=0;
    rend{i}.data = cell(size(dat,1),1);
    if issparse(rend{i}.ren)
        % Assume that images have been DCT compressed
        % - the SPM99 distribution was originally too big.
        d = size(rend{i}.ren);
        B1 = spm_dctmtx(d(1),d(1));
        B2 = spm_dctmtx(d(2),d(2));
        rend{i}.ren = B1*rend{i}.ren*B2';
        % the depths did not compress so well with
        % a straight DCT - therefore it was modified slightly
        rend{i}.dep = exp(B1*rend{i}.dep*B2')-1;
    end
    rend{i}.ren(rend{i}.ren>=1) = 1;
    rend{i}.ren(rend{i}.ren<=0) = 0;
    if showbar, spm_progress_bar('Set', i); end
end
if showbar, spm_progress_bar('Clear'); end

if showbar, spm_progress_bar('Init', length(dat)*length(rend),...
            'Making pictures', 'Number completed'); end

mx = zeros(length(dat),1)+eps;
mn = zeros(length(dat),1);

for j=1:length(dat)
    
    XYZ = dat(j).XYZ;
    t   = dat(j).t;
    dim = dat(j).dim;
    mat = dat(j).mat;
    
    i = 1; % there is only 1 rendering
    
    % transform from Talairach space to space of the rendered image
    %------------------------------------------------------------------
    M1  = rend{i}.M*mat;
    zm  = sum(M1(1:2,1:3).^2,2).^(-1/2);
    M2  = diag([zm' 1 1]);
    M  = M2*M1;
    cor = [1 1 1 ; dim(1) 1 1 ; 1 dim(2) 1; dim(1) dim(2) 1 ;
        1 1 dim(3) ; dim(1) 1 dim(3) ; 1 dim(2) dim(3); dim(1) dim(2) dim(3)]';
    tcor= M(1:3,1:3)*cor + M(1:3,4)*ones(1,8);
    off = min(tcor(1:2,:)');
    M2  = spm_matrix(-off+1)*M2;
    M  = M2*M1;
    xyz = (M(1:3,1:3)*XYZ + M(1:3,4)*ones(1,size(XYZ,2)));
    d2  = ceil(max(xyz(1:2,:)'));
    
    % Calculate 'depth' of values
    %------------------------------------------------------------------
    if ~isempty(d2)
        dep = spm_slice_vol(rend{i}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
        z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));
        
        if ~isfinite(brt), msk = find(xyz(3,:) < (z1+20) & xyz(3,:) > (z1-5));
        else msk = find(xyz(3,:) < (z1+60) & xyz(3,:) > (z1-5)); end
    else
        msk = [];
    end
    
    if ~isempty(msk)
        
        % Generate an image of the integral of the blob values.
        %--------------------------------------------------------------
        xyz = xyz(:,msk);
        if ~isfinite(brt), t0  = t(msk);
        else
            dst = xyz(3,:) - z1(msk);
            dst = max(dst,0);
            t0  = t(msk).*exp((log(0.5)/10)*dst)';
        end
        X0  = full(sparse(round(xyz(1,:)), round(xyz(2,:)), t0, d2(1), d2(2)));
        hld = 1; if ~isfinite(brt), hld = 0; end
        X   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(rend{i}.dep),hld);
        msk = find(X<0);
        X(msk) = 0;
    else
        X = zeros(size(rend{i}.dep));
    end
    
    % Brighten the blobs
    %------------------------------------------------------------------
    if isfinite(brt), X = X.^brt; end
    
    mx(j) = max([mx(j) max(max(X))]);
    mn(j) = min([mn(j) min(min(X))]);
    
    rend{i}.data{j} = X;
    
    if showbar, spm_progress_bar('Set', i+(j-1)*length(rend)); end
    
end

mxmx = max(mx);
mnmn = min(mn);

if showbar, spm_progress_bar('Clear'); end


%-Display
%==========================================================================
if ~spm('CmdLine')
    Fgraph = spm_figure('GetWin', xSPM.title);
    spm_results_ui('Clear',Fgraph);
    
    % if showbar, hght = 2/3; else hght = 0.5; end
    hght = 1;
    % subplot('Position',[0, 0, 1, hght]);
    ax = axes('Parent',Fgraph,...
        'units','normalized',...
        'Position',[0, 0, 1, hght],...
        'Visible','off');
    image(0,'Parent',ax);
    set(ax,'YTick',[],'XTick',[]);
end

rgb = cell(1,length(rend));
% make a parsimanious grid for ploting
ncol = ceil(sqrt(numel(dat)/2));
nrow = 2*ncol;
while (nrow-1)*ncol >= numel(dat); nrow = nrow - 1; end
hght = 0.95;

if ~isfinite(brt)
    % Old style split colourmap display - actually the one that is modified
    % compared to spm_render.m
    %----------------------------------------------------------------------
    split = [gray(64); hot(64)];
    if ~spm('CmdLine'), colormap(split); end
    for i=1:length(dat)
        ren = rend{1}.ren;
        X   = (rend{1}.data{i}-mnmn)/(mxmx-mnmn);
        msk = find(X);
        ren(msk) = X(msk)+(1+1.51/64);
        if ~spm('CmdLine')
            ax = axes('Parent',Fgraph,...
                'units','normalized',...
                'Position',[rem(i-1,ncol)*1/ncol, 1-floor((i-1)/ncol)*hght/nrow-hght/nrow , 1/ncol, hght/nrow],...
                'Visible','off');
            image(ren*64,'Parent',ax);
            set(ax,'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatioMode','auto',...
                'YTick',[],'XTick',[],...
                'XDir','normal','YDir','normal');
            
            % Add subject number
            text('Parent',ax,...
                'Color', [1 1 1],...
                'FontUnits', 'normalized',...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','left',...
                'Position', [1 1],...
                'String', sprintf(' %02.0f', i));
        end
    end

 
else
    % Combine the brain surface renderings with the blobs, and display
    % using 24 bit colour
    %----------------------------------------------------------------------
    for i=1:numel(dat)
        ren = rend{1}.ren;
        X = cell(3,1);
        X{1} = rend{1}.data{i}/(mxmx-mnmn)-mnmn;
        for j=2:3
            X{j}=zeros(size(X{1}));
        end

        rgb{1} = zeros([size(ren) 3]);
        tmp = ren.*max(1-X{1}-X{2}-X{3},0);
        for k = 1:3
            rgb{1}(:,:,k) = tmp + X{1}*col(1,k) + X{2}*col(2,k) +X{3}*col(3,k);
        end
        
        if ~spm('CmdLine')
            ax = axes('Parent',Fgraph,...
                'units','normalized',...
                'Position',[rem(i-1,ncol)*1/ncol, 1-floor((i-1)/ncol)*hght/nrow-hght/nrow , 1/ncol, hght/nrow],...
                'nextplot','add', ...
                'Visible','off');
            image(rgb{1},'Parent',ax);
            set(ax,'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatioMode','auto',...
                'YTick',[],'XTick',[],...
                'XDir','normal','YDir','normal');
            
            % Add subject number
            text('Parent',ax,...
                'Color', [1 1 1],...
                'FontUnits', 'normalized',...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','left',...
                'Position', [1 1],...
                'String', sprintf(' %02.0f', i));
        end
        
        % rgb{i} = flipud(rgb{i}); % changed by FM
    end
end

spm('Pointer','Arrow');

if nargout, varargout = { rgb }; end


%==========================================================================
% function surf_rend(dat,rend,col)
%==========================================================================
function surf_rend(dat,rend,col)

%-Setup figure and axis
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);

ax0 = axes(...
    'Tag',      'SPMMeshRenderBackground',...
    'Parent',   Fgraph,...
    'Units',    'normalized',...
    'Color',    [1 1 1],...
    'XTick',    [],...
    'YTick',    [],...
    'Position', [-0.05, -0.05, 1.05, 0.555]);

ax = axes(...
    'Parent',   Fgraph,...
    'Units',    'normalized',...
    'Position', [0.05, 0.05, 0.9, 0.4],...
    'Visible',  'off');

H = spm_mesh_render('Disp',rend,struct('parent',ax));
spm_mesh_render('Overlay',H,dat,col);
camlight(H.light);

try
    setAllowAxesRotate(H.rotate3d, setxor(findobj(Fgraph,'Type','axes'),ax), false);
end
    
%-Register with MIP
%--------------------------------------------------------------------------
try % meaningless when called outside spm_results_ui
    hReg = spm_XYZreg('FindReg',spm_figure('GetWin','Interactive'));
    xyz  = spm_XYZreg('GetCoords',hReg);
    hs   = mydispcursor('Create',ax,dat.mat,xyz);
    spm_XYZreg('Add2Reg',hReg,hs,@mydispcursor);
end
  

%==========================================================================
function varargout = mydispcursor(varargin)

switch lower(varargin{1})
    %======================================================================
    case 'create'
    %======================================================================
    % hMe = mydispcursor('Create',ax,M,xyz)
    ax  = varargin{2};
    M   = varargin{3};
    xyz = varargin{4};
    
    [X,Y,Z] = sphere;
    vx = sqrt(sum(M(1:3,1:3).^2));
    X = X*vx(1) + xyz(1);
    Y = Y*vx(2) + xyz(2);
    Z = Z*vx(3) + xyz(3);
    hold(ax,'on');
    hs = surf(X,Y,Z,'parent',ax,...
        'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting', 'phong');
    set(hs,'UserData',xyz);
    
    varargout = {hs};
    
    %======================================================================
    case 'setcoords'    % Set co-ordinates
    %======================================================================
    % [xyz,d] = mydispcursor('SetCoords',xyz,hMe,hC)
    hMe  = varargin{3};
    pxyz = get(hMe,'UserData');
    xyz  = varargin{2};
    
    set(hMe,'XData',get(hMe,'XData') - pxyz(1) + xyz(1));
    set(hMe,'YData',get(hMe,'YData') - pxyz(2) + xyz(2));
    set(hMe,'ZData',get(hMe,'ZData') - pxyz(3) + xyz(3));
    set(hMe,'UserData',xyz);
    
    varargout = {xyz,[]};
    
    %======================================================================
    otherwise
    %======================================================================
    error('Unknown action string')

end
