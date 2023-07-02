function gui_lc (action)
%
% This file is used by FASTICAG

% This file holds the callbacks for load-dialog

% 24.8.1998
% Hugo Gävert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables

% Handle to the window
global hf_FastICA_Load;

% Handles to some of the controls in window
global he_FastICA_file;

% Needed handles from the main figure
global ht_FastICA_mixedStatus;

% Needed handles from the advOpt figure
global hb_FastICA_initGuess;
global ht_FastICA_initGuess;
global hpm_FastICA_initState;

% The needed main variables
global g_FastICA_mixedsig;
global g_FastICA_mixedmean;

% What is the load type of load dialog
global g_FastICA_loadType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should not take long...
watchonInFigure = watchon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Load'
 
 varName = get(he_FastICA_file, 'String');      % The name of the variable to be loaded
 command=['evalin(''base'',''assignin(''''caller'''',''''data'''',' varName ')'')'];
 eval(command,'fprintf(''Variable not found in MATLAB workspace, data not loaded!\n'');data=[];');                          % Variable is copyed to 'data'

 if isempty(data)  % if there was no name given...
   watchoff (watchonInFigure);
   %break;
   return;
 end
 switch g_FastICA_loadType
  case  'data'                                % New data
   g_FastICA_mixedsig = data;
   set(ht_FastICA_mixedStatus, 'String', '');
   g_FastICA_mixedmean = [];                 % New data - so that means ...
   gui_cb NewData;                             
   
  case 'guess'                                % New initial guess
   set(hb_FastICA_initGuess, 'UserData', data);     % Since we loaded new initial
   set(ht_FastICA_initGuess, 'String', 'Loaded');   % guess, we wan't to use it too
   set(hpm_FastICA_initState, 'Value', 2);          % ... set initState to 'guess'
 end

 close(hf_FastICA_Load);                       % close the dialog

 % Use break to 'jump' over the watchoff statement at the end
 %break;
 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Cancel'
 
 close(hf_FastICA_Load);                       % do nothing just exit

 % Use break to 'jump' over the watchoff statement at the end
 %break;
 return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'Help'

 % Which help do we show?
 switch g_FastICA_loadType
  case 'data'
   gui_help('gui_lc_data');
  case 'guess'
   gui_help('gui_lc_guess');
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    % switch

watchoff (watchonInFigure);
