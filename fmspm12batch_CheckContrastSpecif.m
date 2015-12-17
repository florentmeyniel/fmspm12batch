function fmspm12batch_CheckContrastSpecif(modelname)
% Automatically get the contrast specification, with names from the
% multicond. This is for sanity checking.
% 
% Usage: fmspm12batch_CheckContrastSpecif(modelname)

% Get parameter definition
eval(sprintf('%s_paramfile', modelname))

for iSub = 1:numel(sublist)
    fprintf('\n\nSUBJECT %d', sublist(iSub))
    
    % get contrast file
    fcon = sprintf('%s/subj%02.0f/MultiCond/%s_contrast.mat', datadir, sublist(iSub), modelname);
    dcon = load(fcon);
    
    % get multicond file
    fmcond = sprintf('%s/subj%02.0f/MultiCond/%s_multicond_session%d.mat', datadir, sublist(iSub), modelname, 1);
    dmcond = load(fmcond);
    
    % get number of sessions
    flist = dir(sprintf('%s/subj%02.0f/MultiCond/%s_multicond_session*.mat', datadir, sublist(iSub), modelname));
    nSess = numel(flist);
    fprintf('\nNUMBER OF SESSION IDENTIFIED: %d', nSess)
    
    % get list of regressors from the 1st session
    nRegOns = numel(dmcond.names);
    iCount = 0;
    RegName = {};
    for iReg = 1:nRegOns
        iCount = iCount + 1;
        RegName{iCount} = dmcond.names{iReg};
        if ~isfield(dmcond.pmod(iReg), 'name') || isempty(dmcond.pmod(iReg).name)
            nRegPmod = 0;
        else
            nRegPmod = numel(dmcond.pmod(iReg).name);
        end
        for iRegPmod = 1:nRegPmod
            iCount = iCount + 1;
            RegName{iCount} = [dmcond.names{iReg}, ' * ', dmcond.pmod(iReg).name{iRegPmod}];
        end
    end
    
    % complete with number of sessions, movement parameters & constants
    fullRegName = {};
    for iSess = 1:nSess
        if strcmp(bases_functions, 'hrf')
            fullRegName = {fullRegName{:}, RegName{:}, ...
                'mov par', 'mov par', 'mov par', 'mov par', 'mov par', 'mov par'};
        elseif strcmp(bases_functions, 'hrf+deriv')
            for k = 1:numel(RegName)
                fullRegName = {fullRegName{:}, RegName{k}, 'deriv'};
            end
            fullRegName = {fullRegName{:}, ...
                'mov par', 'mov par', 'mov par', 'mov par', 'mov par', 'mov par'};
        end
    end
    for iSess = 1:nSess
        fullRegName = {fullRegName{:}, 'cst'};
    end
        
    % Get matches in contrast definition
    for iCon = 1:numel(dcon.names)
        % Get Contrast Name
        fprintf('\n%s', dcon.names{iCon})
        
        % Get matches
        ind = find(dcon.values{iCon});
        for k = 1:numel(ind)
            if dcon.values{iCon}(ind(k)) > 0
                fprintf('\n\t(+) %s', fullRegName{ind(k)})
            else
                fprintf('\n\t(-) %s', fullRegName{ind(k)})
            end
        end
        
        % prepare a warning message if the contrast is biased.
        warningmess = '';
        if sum(abs(dcon.values{iCon})) ~= 1
            warningmess = '!!!';
        else
            if ~all(dcon.values{iCon}>0) && ~all(dcon.values{iCon}<0)
            else
                if sum(dcon.values{iCon}) ~= 0
                    warningmess = '!!!';
                end
            end
        end
           
        fprintf('\n\t     sum(c):%d sum(abs(c)):%d %s', ...
            sum(dcon.values{iCon}), sum(abs(dcon.values{iCon})), warningmess)
    end
end