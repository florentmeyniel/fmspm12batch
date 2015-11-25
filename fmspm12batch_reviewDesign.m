function fmspm12batch_reviewDesign(spmname)
% Function to review the design specified in a SPM.mat:
%   * Show regressors in the CheckOrth_fn style
%   * Compute mean statistical power
%   * plot the mean design orthogonality
%
% NB: assume that session are identical.
%
% Usage: fmspm12batch_reviewDesign(spmname)
%   spmname: name of the SPM.mat to load

% load SPM.mat to review
load(spmname)

% Review regressors visually
% =========================================================================
try
    CheckOrth_fn(SPM, 1, 1, 'w')
catch
    CheckOrth_fn(SPM, 1, 1, 'o')
end

% Power analysis
% =========================================================================

% 1- Compute power
% ~~~~~~~~~~~~~~~~
% NB: assume that all sessions have the same regressors

nSess = numel(SPM.Sess);

% Get filtered & whitened design matrix
% NB: the formula with inv(X'X) works for a filtered & whitened design
% matrix, see SPM book, Chap 15 by Henson on 'Efficient Experimental
% Design for fMRI', eq. 15.2 to 15.4
X = SPM.xX.xKXs.X(:,1:end-nSess);

% turn the 0 padding a NaN to remove it from computations
XNaN = ones(size(X));
XNaN(X == 0) = NaN;
X = X .* XNaN;

% de-mean design matrix to compute co-variance matrix easily
X = (X - ones(size(X,1), 1)*nanmean(X));

% replace the NaN padding back by a 0 padding
X(isnan(X)) = 0;

% compute the inverse of co-variance matrix
V = inv((X'*X));
nReg = size(X, 2);
E = zeros(1, nReg);
for k = 1:nReg
    c = zeros(nReg, 1);
    c(k) = 1;
    E(k) = 1/(c'*V*c);
end

% sum over sessions
reshape(E, [nSess, nReg/nSess]);
statpow = sum(reshape(E, [nSess, nReg/nSess]));

% remove movement parameters
statpow = statpow(1:end-6);

% Get list of regressors
% ~~~~~~~~~~~~~~~~~~~~~~
% get list of regressors from the 1st session
nRegOns = numel(SPM.Sess(1).U);
iCount = 0;
RegName = {};
isOns = [];
isPmod = [];
for iReg = 1:nRegOns
    iCount = iCount + 1;
    
    % get type
    isOns(iCount) = 1;
    isPmod(iCount) = 0;
    
    RegName{iCount} = SPM.Sess(1).U(iReg).name{1};
    nRegPmod = numel(SPM.Sess(1).U(iReg).name)-1;
    for iRegPmod = 1:nRegPmod
        iCount = iCount + 1;
        
        % get type
        isOns(iCount)   = 0;
        isPmod(iCount) = 1;
        
        RegName{iCount} = SPM.Sess(1).U(iReg).name{1+iRegPmod};
    end
end
isOns = logical(isOns);
isPmod = logical(isPmod);

if SPM.xBF.order == 1
    ind = ones(1, numel(RegName));
elseif SPM.xBF.order == 2
    ind = repmat([1 0], [1, numel(RegName)]);
elseif SPM.xBF.order == 3
    ind = repmat([1 0 0], [1, numel(RegName)]);
end
ind = logical(ind);

% 3- plot results
% ~~~~~~~~~~~~~~~~~~~~~~
figure;
set(gcf, 'Name', 'Statistical Power')
set(gcf, 'Color', [1 1 1])
subplot(2, 1, 1)
tmp = statpow(ind);
barh(tmp(isOns))
set(gca, 'YTick', 1:sum(isOns), 'YTickLabel', {RegName{isOns}})

subplot(2, 1, 2)
tmp = statpow(ind);
h = barh(tmp(isPmod));
set(gca, 'YTick', 1:sum(isPmod), 'YTickLabel', {RegName{isPmod}})

% Check design othogonality
% =========================================================================

% Get computation from spm_DesRep.m
[nScan,nPar] = size(SPM.xX.xKXs.X);
% taken from spm_DesRep, l. 679 to 694
tmp = sqrt(sum(SPM.xX.xKXs.X.^2));
O   = SPM.xX.xKXs.X'*SPM.xX.xKXs.X./kron(tmp',tmp);
tmp = sum(SPM.xX.xKXs.X);
tmp     = abs(tmp)<eps*1e5;

% taken from spm_DesRep, l. 799
tmp = 1-abs(O); 
tmp(logical(tril(ones(nPar),-1))) = 1;

Ortho = tmp;

% remove constant
Ortho = Ortho(1:end-nSess, 1:end-nSess);

% Average over session
% ~~~~~~~~~~~~~~~~~~~~
nCol = size(Ortho,2);
nReg_sess = nCol/nSess;
meanOrtho = zeros(nReg_sess);
for iSess = 1:nSess;
    meanOrtho = meanOrtho + (1/nSess)*...
        Ortho(nReg_sess*(iSess-1)+1:nReg_sess*(iSess), ...
        nReg_sess*(iSess-1)+1:nReg_sess*(iSess));
end

% Plot result
% ~~~~~~~~~~~
figure;
set(gcf, 'Name', 'Design Orthogonality')
set(gcf, 'Color', [1 1 1])
imagesc(meanOrtho)
axis square
colormap('bone')
set(gca, 'YTick', 1:SPM.xBF.order:SPM.xBF.order*numel(RegName), ...
    'YTickLabel', {RegName{:}})
