function fmspm12batch_replotRealign(dirname, spm_dir)
% Function to load the files with realignment parameters, and plot them in
% the 'SPM style'.
% 
% fpmsp12batch_replotRealign(dirname, spm_dir)
%   dirname: path to folder where the realignment files can be found
%   spm_dir: root directory of spm.

% initialize spm
% =========================================================================

close all
addpath(spm_dir)

% get data
% =========================================================================

% list file name (ascii files with realignement parameters)
filename = cellstr(spm_select('FPList', dirname, '^rp_aepi_.*\.txt'));

% load data from ascii files
Params = [];
for k = 1:numel(filename)
    tmp = load(filename{k});
    Params = [Params; tmp];
end

% plot data
% =========================================================================

% (adapted from plot_parameters, in the spm_realign.m)
fg = spm_figure('Create','Graphics');
ax = axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'Visible','off');
set(get(ax,'Title'),'String','Image realignment',...
    'FontSize',16,'FontWeight','Bold','Visible','on');

ax = axes('Position',[0.1 0.35 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Params(:,1:3),'Parent',ax)
s  = {'x translation','y translation','z translation'};
%text([2 2 2], Params(2, 1:3), s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','mm');


ax = axes('Position',[0.1 0.05 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Params(:,4:6)*180/pi,'Parent',ax)
s  = {'pitch','roll','yaw'};
%text([2 2 2], Params(2, 4:6)*180/pi, s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','degrees');

