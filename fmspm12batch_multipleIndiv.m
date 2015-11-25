% Wrapper function to display on the same figure the individual statistics.
% The script prompt for a display mode and run the corresponding function.
% Available modes are:
%   MIP
%   Slices
%   Surface


% Initialize SPM
try ver = spm('Version'); 
catch
    error('the SPM toolbox is not in the Matlab''s path')
end
spm_results_ui('clear')
spm_results_ui('close')

% Select a view
ViewType = spm_input('Display multiple subjects as', '+1', ...
    'MIP|Slices|Surface', 1:3, 1);

% Run the scrip corresponding to the requested view
if ViewType == 1
    fmspm12batch_multipleMIP
elseif ViewType == 2
    fmspm12batch_multipleSlices
elseif ViewType == 3
    fmspm12batch_multipleRender
end