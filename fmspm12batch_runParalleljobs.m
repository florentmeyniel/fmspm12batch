function fmspm12batch_runParalleljobs(matfile, spm_path, nJobMax, MatlabCmd)
% Function to launch batches in parallel, from separate matlab workers
%
% fmspm12batch_runParalleljobs(matfile, spm_path, nJobMax, MatlabCmd))
%   matfile: structure with fullpath of mat file describing the SPM job
%   spm_bath: path of the SPM toolbox
%   nJobMax: maximal number of independent workers that are allowed at any
%       given moment.
%   MatlabCmd: command to invoke Matlab from a terminal
%
% Note: if the number of Matlab open exceeds the number of CPU, the parallelization
% becomes stupidely sub-optimal, with only one CPU used...

% CHECK THAT ENOUGH CPUs ARE AVAILABLE
[~, w] = unix('top -n 1 | grep MATLAB');
Nopen = numel(strfind(w, 'MATLAB'));
[~, w] = unix('nproc');
Ncpu = str2num(w);
if Nopen > Ncpu
    fprintf('\n ***************************************************')
    fprintf('\n 	              WARNING                      ')
    fprintf('\n     There are more Matlab (%d) than cores (%d)     ', Nopen, Ncpu)
    fprintf('\n   The parallelization is likely to be inefficient  ')
    fprintf('\n ***************************************************')
    fprintf('\n')
end

% Get the fmspm12batch patch
currpath = mfilename('fullpath');
ind = strfind(currpath, '/');
toolboxpath = currpath(1:ind(end));

% Initialize list of jobs
nJob = numel(matfile);
JobList = 1:nJob;
RunningJobs = [];

% get job type
JobType = matfile{1};
ind = strfind(JobType, '_');
JobType = JobType(ind(end)+1:end-4);

% Launch jobs in parallel, given the maximum number of job allowed in
% parallel
timestamp = sprintf('%d-%d-%d_%d-%d-%1.0f_%5.0f', clock, rand*1e5);
while ~isempty(JobList)
    
    % LAUNCH A BATCH IF ALLOWED
    if length(RunningJobs) < nJobMax
        % Get job number
        iJob = JobList(1);
        RunningJobs = [RunningJobs, iJob];
        
        % udpate the list of remaining jobs
        if numel(JobList) > 1
            JobList = JobList(2:end);
        else
            JobList = [];
        end
        
        if strmatch(JobType, 'TopupCorrection')
            cmdstr = [sprintf('%s -nosplash -nodesktop -r "', MatlabCmd), ...           % open Matlab
                sprintf('tic; '), ...                                                   % tic
                sprintf('addpath %s ; ', toolboxpath), ...                              % go into correct directory
                sprintf(['try fmspm12batch_AddTopupCorrection_job1sub(''%s'');' ...
                         ' catch ME; disp(getReport(ME,''extended'')); end ; '], ...    % execute function of job
                matfile{iJob}) ...
                sprintf('mkdir(''over_%s_%s_%d'') ; ', timestamp, JobType, iJob), ...   % make a temporary directory to signal job end
                sprintf('toc; '), ...                                                   % toc
                'quit' ...                                                              % close Matlab
                '" '...
                sprintf('> log_file_%s_%d.txt ', JobType, iJob), ...
                '&'];
        else
            
            % make command to exectute this job
            cmdstr = [sprintf('%s -nosplash -nodesktop -r "', MatlabCmd), ...           % open Matlab
                sprintf('tic; '), ...                                                   % tic
                sprintf('addpath %s ; ', toolboxpath), ...                              % go into correct directory
                sprintf(['try fmspm12batch_run1job(''%s'', ''%s'');'...
                         ' catch ME; disp(getReport(ME,''extended'')); end ; '], ...    % execute function of job
                matfile{iJob}, spm_path) ...
                sprintf('mkdir(''over_%s_%s_%d'') ; ', timestamp, JobType, iJob), ...   % make a temporary directory to signal job end
                sprintf('toc; '), ...                                                   % toc
                'quit' ...                                                              % close Matlab
                '" '...
                sprintf('> log_file_%s_%d.txt ', JobType, iJob), ...
                '&'];
        end
        
        % execute command
        system(cmdstr);
        
        % print job name in command line window
        fprintf('\n Job launched: %s', matfile{iJob})
    end
    
    % UPDATE LIST OF BATCHES CURRENTLY RUNNING
    for iJob = RunningJobs
        % check if job has returned (with the 'probe directory')
        if exist(sprintf('over_%s_%s_%d', timestamp, JobType, iJob), 'dir')
            % remove the probe
            rmdir(sprintf('over_%s_%s_%d', timestamp, JobType, iJob))
            
            % Remove job from list of running jobs
            RunningJobs = setdiff(RunningJobs, iJob);
            fprintf('\n Job completed: %s', matfile{iJob})
        end
    end
    
    % PAUSE THE FUNCTION
    % after nJobMax jobs have been launched
    % to avoid running the while loop frenetically...
    if max(RunningJobs) >= nJobMax
        pause(1)
    end
end

fprintf('\n Check that all jobs are done to return')
% When all jobs have been launched, wait until they are all terminated to
% return (and print in the command window that jobs have returned)
while ~isempty(RunningJobs)
    
    % UPDATE LIST OF BATCHES CURRENTLY RUNNING
    for iJob = RunningJobs
        % check if job has returned (with the 'probe directory')
        if exist(sprintf('over_%s_%s_%d', timestamp, JobType, iJob), 'dir')
            % remove the probe
            rmdir(sprintf('over_%s_%s_%d', timestamp, JobType, iJob))
            
            % Remove job from list of running jobs
            RunningJobs = setdiff(RunningJobs, iJob);
            fprintf('\n Job completed: %s', matfile{iJob})
        end
    end
    
    % PAUSE THE FUNCTION
    % to avoid running the while loop frenetically...
    pause(1)
end

fprintf('\n')
 
