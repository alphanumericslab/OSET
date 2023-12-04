function gui_sc (action)
%
% This file is used by FASTICAG

% This file holds the callbacks for save-dialog

% 5.6.1999
% Hugo Gävert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the window
global hf_FastICA_Save;

% Handles to some of the controls in window
global he_FastICA_suffix;

% The needed main variables
global g_FastICA_ica_sig;
global g_FastICA_ica_A;
global g_FastICA_ica_W;
global g_FastICA_white_sig;
global g_FastICA_white_wm;
global g_FastICA_white_dwm;
global g_FastICA_pca_E;
global g_FastICA_pca_D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should not take long...
watchonInFigure = watchon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Save'
 
 suffix = deblank(get(he_FastICA_suffix, 'String')); % The suffix for the variables

 fprintf('Saving results in variables in Matlab workspace.\n');
 assignin('base',['IC' suffix],g_FastICA_ica_sig);
 assignin('base',['A' suffix],g_FastICA_ica_A);
 assignin('base',['W' suffix],g_FastICA_ica_W);
 assignin('base',['whitesig' suffix],g_FastICA_white_sig);
 assignin('base',['whiteningMatrix' suffix],g_FastICA_white_wm);
 assignin('base',['dewhiteningMatrix' suffix],g_FastICA_white_dwm);
 assignin('base',['E' suffix],g_FastICA_pca_E);
 assignin('base',['D' suffix],g_FastICA_pca_D);

 close(hf_FastICA_Save);                  % close the dialog
 
 % Use break to 'jump' over the watchoff statement at the end
 break;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Cancel'
 
 close(hf_FastICA_Save);                       % do nothing just exit

 % Use break to 'jump' over the watchoff statement at the end
 break;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Help'

 gui_help('gui_sc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    % switch

watchoff (watchonInFigure);
