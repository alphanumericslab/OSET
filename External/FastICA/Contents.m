% FastICA for Matlab 5.x
% Version 2.1, January 15 2001
% Copyright (c) Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen.
%
% Type fasticag to launch the graphical user interface
%
% Please refer to your Matlab documentation on how to add FastICA to your
% Matlab search path. (One place to start is the path-command)
%
% FastICA programs:
%   fasticag  - Graphical user interface for FastICA
%   fastica   - command line version of FastICA
%
% Separate functions used by FastICA programs.
%   fpica     - main algorithm for calculating ICA
%   whitenv   - function for whitening data
%   pcamat    - calculates the PCA for data
%   remmean   - function for removing mean
%
%   gui_cb    - needed by fasticag
%   gui_adv   - needed by fasticag
%   gui_advc  - needed by fasticag
%   gui_l     - needed by fasticag
%   gui_lc    - needed by fasticag
%   gui_s     - needed by fasticag
%   gui_sc    - needed by fasticag
%   gui_cg    - needed by fasticag
%   gui_help  - needed by fasticag
%
%   icaplot   - for plotting the signals
%               (also used by fastica and fasticag)
%
% Misc.
%   demosig   - generates some test signals
%
% Deprecated
%   dispsig   - plots the data vectors
%               replaced by icaplot