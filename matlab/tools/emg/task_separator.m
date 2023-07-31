function sep_file_names = task_separator(filename)
% 
% sep_file_names = task_separator(filename)
% 
% *************************************************************************
% * Task based separation of TDSKA/EEG data to avoid memory overuse       *                                      *
% *************************************************************************
% 
% Usage:    sep_file_names = task_separator(filename)
% inputs:
%           'filename': original whole-data file name as a string
% outputs:
%           'sep_file_names': cell array containing separated files' names 
%                             for later use in either 'test_TDSKAEEG_ERPanalysis_main.m'
%                             or 'test_TDSKAEEG_TFanalysis_main.m'
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
    error('***wrong number of input arguments. Refer to manual for details***')
else
    if ~ischar(filename)
        error('***file name has to be a string***')
    end
end

%% initializing & loading data
fprintf('\n Starting the Separation Process...\n');
fprintf('1.5 Loading Data...\n');
data = load(filename);
data = data.data; m = length(data);

%% separating data, saving and returning file names
fprintf('2.5 Separating Data and Saving New Files...\n');
sep_file_names = cell(1, m);
for i=1:m
    task = data{i};
    [task_name, ~] = strtok(task.dataname, '.');
    new_filename = strcat('data_', task_name);
    save(new_filename, 'task');
    sep_file_names{i} = new_filename;
end
fprintf('3.5 Files Saved..!\n');
save(strcat(filename, '-filenames'), 'sep_file_names')
fprintf('4.5 New Files Names Recorded..!\n');

%% clean-up memory
clearvars -except sep_file_names
fprintf('5.5 Memory Cleaned Up..!\n');
fprintf('Done..!\n');

end
