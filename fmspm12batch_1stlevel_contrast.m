function consess = fmspm12batch_1stlevel_contrast(datadir, iSub, modelname)        
% Function to specify contrasts. 
%
% consess = fmspm12batch_1stlevel_contrast(datadir, iSub, modelname)
%  	Outputs:
% 		* consess: matlabbatch like structure with contrast specification
% 	Inputs:
% 		* iSub: subject number, to find the subject folder subjXX
% 		* dataloc: structure to locate data
% 		* modelname: name of the model.
% 
% The input specify where the .mat with the contrast specification should be found:
% %s/subj%02.0f/MultiCond/%s.mat ; datadir, iSub, modelname
%


% load the contrast file
fname = sprintf('%s/subj%02.0f/MultiCond/%s_contrast.mat', datadir, iSub, modelname);
dat = load(fname);
        

% Specify T-contrast
nCon = numel(dat.names);
consess = {};
for iCon = 1:nCon
    if strcmp(dat.type{iCon}, 'T')
        consess{iCon}.tcon.name     = dat.names{iCon};
        consess{iCon}.tcon.weights  = dat.values{iCon};
        consess{iCon}.tcon.sessrep  = 'none';
    elseif strcmp(dat.type{iCon}, 'F')
        consess{iCon}.fcon.name     = dat.names{iCon};
        consess{iCon}.fcon.weights  = dat.values{iCon};
        consess{iCon}.fcon.sessrep  = 'none';
    else
        error('Contrast type %d not recognized', dat.type{iCon})
    end
end
