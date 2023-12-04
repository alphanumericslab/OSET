function [Tr_sig, loc] = sig_trend(sig)
% 
% [Tr_sig, loc] = sig_trend(sig)
% 
% *************************************************************************
% * Calculating the trend of a signal using local minima                  *
% *************************************************************************
% 
% Usage:    [Tr_sig, loc] = sig_trend(sig)
% inputs:
%           'sig': input raw signal
% outputs:
%           'Tr_sig': trend vector
%           'loc': returns the locations required for plotting the trend
%
% This program is provided by ESMAEIL SERAJ. Please make sure to cite BOTH 
% the original studies and the User Manual to help others find these items.
% 
% Authors:
% 			Esmaeil Seraj, Karthiga Mahalingam
% Websites:
%			https://github.com/EsiSeraj/ERP_Connectivity_EMG_Analysis
% 			http://oset.ir/category.php?dir=Tools
% 
%    Copyright (C) <2018>  <ESMAEIL SERAJ (eseraj3@gatech.edu)>
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program. If not, see <http://www.gnu.org/licenses/> or
%    write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth
%    Floor, Boston, MA  02110-1301, USA.
% 

%% Checking inputs and assigning default values
if nargin ~= 1
    error('***wrong number of input arguments. Refer to Manual for details***')
else
    if isscalar(sig)
        error('***input has to be a double vector or matrix***')
    end
end
    
%% trend extraction
[maxima, ~] = findpeaks(sig);
[~, loc_minima] = findpeaks(-sig);
minima = sig(loc_minima);
m = min(length(maxima), length(minima));
maxima = maxima(1:m);
minima = minima(1:m);
loc_minima = loc_minima(1:m);
Tr_sig = (maxima + minima)/2;
loc = loc_minima;

end