% Script to move and rename physiological file.
% Before using this script:
%   - put your physiological files in root1/subj01, root/subj02, etc.
%   - rename the file names so that they are meaningful (replace the
%     horrible native name by "sess1", "sess2" for instance).
%   - put you fMRI files in root2/subj01/fMIR, root2/subj02/fMIR, etc.

% This script assume that the fMRI files have names like
% epi_sess1_*.nii, epi_sess2_*.nii, etc. and that the physiological files
% have names like *sess1_Info.log, *sess2_Info.log, etc. Info.log is a
% particular suffix, a list can be provided.
% All files containing the specified suffix and session number will be
% rename and more to the folder containing the fMRI files.

% If files do not include sess1, sess2, etc. and are organized differently,
% tweaking this code will be required.

% add SPM to the path
addpath /volatile/meyniel/toolbox/matlab/spm12/

% Where the subject-specific folder with physiological files are.
source_dir = '/neurospin/unicog/protocols/IRMf/MeynielMoreno_NACONF_2016/Physio_data/raw_data/';

% Where the subject-specific folders with physiological files should be move.
target_dir = '/neurospin/unicog/protocols/IRMf/MeynielMoreno_NACONF_2016/MRI_data/raw_data';

% List of subject
sub_list = 3:7; 

% list of session to process
sess_list = 1:4;

% list of suffices to process
suffix_list = {'Info.log', 'PULS.log', 'RESP.log'};


% RENAME AND MOVE
for iSub = 1:numel(sub_list)
    for iSess = 1:numel(sess_list)
        for iSf = 1:numel(suffix_list)
            
            % find the file that contains this session number and suffix
            datdir = sprintf('%s/subj%02.0f/FileNameWithSessionNumber/', source_dir, sub_list(iSub));
            fname_phys = spm_select('List', datdir, sprintf('^.*sess%d_%s', sess_list(iSess), suffix_list{iSf}));
            
            if isempty(fname_phys)
                fprintf('\n Subj %d, sess %d, suffix %s not found', ...
                    sub_list(iSub), sess_list(iSess), suffix_list{iSf})
            elseif size(fname_phys,1)>1
                fprintf('\n Subj %d, sess %d, suffix %s, %d files found ', ...
                    sub_list(iSub), sess_list(iSess), suffix_list{iSf}, size(fname_phys,1))
            else
                
                % find the corresponding epi session
                datdir = sprintf('%s/subj%02.0f/fMRI/', target_dir, sub_list(iSub));
                fname_epi = spm_select('List', datdir, sprintf('^epi_sess%d_.*.nii', sess_list(iSess)));
                
                if isempty(fname_epi)
                    fprintf('\n Subj %d, sess %d, nii file not found', ...
                        sub_list(iSub), sess_list(iSess))
                elseif size(fname_epi,1)>1
                    fprintf('\n Subj %d, sess %d, %d nii files found ', ...
                        sub_list(iSub), sess_list(iSess), size(fname_epi,1))
                else
                    
                    % move file to target directory
                    source_fname = sprintf('%s/subj%02.0f/FileNameWithSessionNumber/%s', ...
                        source_dir, sub_list(iSub), fname_phys);
                    target_fname = sprintf('%s/subj%02.0f/fMRI/%s_%s', ...
                        target_dir, sub_list(iSub), fname_epi(1:end-4), suffix_list{iSf});
                    copyfile(source_fname, target_fname)
                end
            end
            
        end
    end
end
