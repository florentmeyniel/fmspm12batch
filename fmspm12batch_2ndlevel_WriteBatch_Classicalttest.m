function [matfile, matlabbatch] = fmspm12batch_2ndlevel_WriteBatch_Classicalttest(dataloc, param, confilelist)
% Function to specify the 2nd-level t-test.
%
% [matfile, matlabbatch] = fmspm12batch_2ndlevel_WriteBatch_Classicalttest(dataloc, param, confilelist)
%  	Outputs:
% 		* matfile: full path of the batch saved.
% 		* matlabbatch: batch structure (that is saved)
% 	Inputs:
% 		* dataloc: structure to locate data
% 		* param: structure with specification parameters.
% 		* confilelist: cell with the contrast files to take into a t-test
%


% Specify design matrix
% =========================================================================

matlabbatch{1}.spm.stats.factorial_design.dir                       = {dataloc.con_dir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = confilelist;
matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;                   % implicit masking
matlabbatch{1}.spm.stats.factorial_design.masking.em                = {param.mask};        % explicit masking
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;                   % no Grand Mean Scaling
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;                   % no normalization

matlabbatch{2}.spm.stats.fmri_est.spmmat(1)                         = {sprintf('%s/SPM.mat', dataloc.con_dir)};
matlabbatch{2}.spm.stats.fmri_est.write_residuals                   = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical                  = 1;                    % classical t-test

matlabbatch{3}.spm.stats.con.spmmat(1)                              = {sprintf('%s/SPM.mat', dataloc.con_dir)};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name                   = sprintf('(+) %s', dataloc.con_name);
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights                = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep                = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name                   = sprintf('(-) %s', dataloc.con_name);
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights                = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep                = 'none';
matlabbatch{3}.spm.stats.con.delete                                 = param.delprevcon;

% Save SPM batch
% =========================================================================

matfile = sprintf('%s/batch_Contrast2ndLevel.mat', dataloc.con_dir);
if exist(matfile) 
    delete(matfile); fprintf('\nremove previous batch\n')
end

fprintf('\n Writing SPM batch at: \n %s\n', matfile)
save(matfile,'matlabbatch');

