function matlabbatch = fmspm12batch_2ndLevelCorr(dataloc, regressor, param)


% Specify the design matrix
matlabbatch{1}.spm.stats.factorial_design.dir                       = {dataloc.con_dir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans            = dataloc.flist; 
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c           = regressor.val(:); 
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname       = regressor.varname; %%
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC         = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint           = 1;
matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em                = {param.mask};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;

% Estimate the design matrix
matlabbatch{2}.spm.stats.fmri_est.spmmat(1)                         = {sprintf('%s/SPM.mat', dataloc.con_dir)};
matlabbatch{2}.spm.stats.fmri_est.write_residuals                   = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical                  = 1;                    % classical t-test

% Report the contrast
matlabbatch{3}.spm.stats.con.spmmat(1)                              = {sprintf('%s/SPM.mat', dataloc.con_dir)};

matlabbatch{3}.spm.stats.con.consess{1}.tcon.name                   = sprintf('(+) %s', dataloc.con_name);
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights                = [0 1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep                = 'none';

matlabbatch{3}.spm.stats.con.consess{2}.tcon.name                   = sprintf('(-) %s', dataloc.con_name);
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights                = [0 -1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep                = 'none';

matlabbatch{3}.spm.stats.con.consess{3}.tcon.name                   = '(+) intercept';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights                = [1 0];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep                = 'none';

matlabbatch{3}.spm.stats.con.consess{4}.tcon.name                   = '(-) intercept';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights                = [-1 0];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep                = 'none';

matlabbatch{3}.spm.stats.con.delete                                 = param.delprevcon;