function fmspm12batch_run1job(jobfile, spm_path)
% Run an SPM job
% syntax: fmspm12batch_run1job(jobfile, spm_path)
% 	jobfile: maybe a matlab cell, or file name (string) of its .mat 
%            (e.g. 'toto.mat' if it is in the current working directory)
% 	spm_path (optional): SPM toolbox directory
% 
% NB: to avoid SPM figures, this function may be run in 
% a terminal after launching matlab without the interface:
% matlab-R2013a -nosplash -nodesktop
% (unfortunately, new versions of Matlab are more restrictive
% on the -nojvm mode, so that the SPM figures make this function
% crash under this mode.)

% check that SPM12 is in the path, or can be set in the path.
if nargin == 1
    if exist('/volatile/meyniel/toolbox/matlab/spm12/')
        addpath /volatile/meyniel/toolbox/matlab/spm12/
    else
        try ver = spm('Version');
        catch
        error('SPM12 is not in the path')
        end
    
    	if ~strcmp(ver(1:5), 'SPM12')
    	    error('the SPM in the path is not SPM12')
    	end
    end
else
    addpath(spm_path)
end

% Initialize spm
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% check that job exist
if iscell(jobfile) 
    fprintf('\n job is specified in a cell\n')
else
    if ischar(jobfile)
        fprintf('\n job is specified in a file:\n%s\n', jobfile)
        if ~exist(jobfile)
            error(sprintf('this jobfile does not exist: \n%d\n', jobfile))
        end
    else
        error('jobfile must be a matlab cell or string')
    end
end

% launch batch
spm_jobman('run', jobfile);
