function toto
% start 2 cmd in matlab, calling a function. 
% The two commands are launched sequentially, and set in the background to 
% run in parallel.
% In addition: print a log of matlab action in the folder where matlab is launched

msg = {'message1', 'message2'};

for imes = 1:length(msg)
    timestamp = sprintf('%d-%d-%d_%d-%d-%1.0f', clock);
    cmdstr = ['matlab-R2013a -nosplash -nodesktop -r "' ...                             % open Matlab
        'cd ~/PostDoc/manip/MarkovConfidencefMRI/script/MRI_script/fmspm12batch ; ' ... % go into correct directory
        sprintf('try ParallelBatchExample_subfn(''%s''); catch; end ; ', msg{imes}) ... % execute function of job
        'quit' ...                                                                      % close Matlab
        '" '...
        sprintf('> log_file_%s_%s.txt ', msg{imes}, timestamp), ...
    	'&'];
    unix(cmdstr)
end
