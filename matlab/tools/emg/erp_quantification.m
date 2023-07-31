function [erd_area, ers_area, quant_erp] = erp_quantification(erp, fs, trigger_time, varargin)
% 
% [erd_area, ers_area, quant_erp] = erp_quantification(erp, fs, trigger_time)
% [erd_area, ers_area, quant_erp] = erp_quantification(erp, fs, trigger_time, ref_per, cof_intv)
% 
% *************************************************************************
% * ERP Area (ERD and ERS events' area) Quantification                    *
% *************************************************************************
% 
% Usage:    [erd_area, ers_area, quant_erp] = erp_quantification(erp, fs, trigger_time)
%           [erd_area, ers_area, quant_erp] = erp_quantification(erp, fs, trigger_time, ref_per, cof_intv)
% inputs:
%           'erp': vector of estimated ERP signal
%           'fs': sampling frequency (Hz)
%           'trigger_time': trigger onset flag (Seconds)
%     (opt) 'ref_per': double vector in [-a, -b] form where -a and -b are 
%                      the edges of reference segment (default: ref_per = 
%                      [-1.3, -0.3] Seconds)
%     (opt) 'cof_intv': confidence interval coefficient (default: cof_intv = 3)
% outputs:
%           'erd_area': ERD events area
%           'ers_area': ERS events area
%           'quant_erp': erp with quantified magnitude
% Note:
%           an empty bracket [] Must be assigned to not-specified values
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
if nargin < 3
    error('***wrong number of input arguments. Refer to Manual for details***')
elseif nargin == 3
    ref_per = [-1.3, -0.3];
    cof_intv = 3;
elseif nargin > 3
    if size(varargin, 2) ~= 2
        error('***an empty bracket [] Must be assigned to not-specified values***')
    else
        if isempty(varargin{1})
            ref_per = [-1.3, -0.3];
        else
            ref_per = varargin{1};
        end
        if isempty(varargin{2})
            cof_intv = 3;
        else
            cof_intv = varargin{2};
        end
    end
end
if (~(isscalar(fs) && isscalar(trigger_time) && isscalar(cof_intv)))
    error('***fs, trigger time and confidence interval coeff. have to be scalars***')
end
if(isscalar(erp) || isscalar(ref_per))
    error('***input ERP signal and the reference period have to be double vectors***')
end

%% calculating the mean and confidence intervals based on reference period
time_vect = round((trigger_time+ref_per(1)))*fs:round((trigger_time+ref_per(2))*fs);
ref_val = mean(erp(time_vect));
intv_tresh = [-cof_intv*sqrt(std(erp(time_vect))), cof_intv*sqrt(std(erp(time_vect)))];

%% quantifying ERP magnitude
quant_erp = (erp/ref_val)*100;

%% separating ERD and ERS events
erp_after = erp(trigger_time*fs:end);
erd_indx = find(erp_after<=(ref_val+intv_tresh(1)));
ers_indx = find(erp_after>=(ref_val+intv_tresh(2)));

%% separating different ERD/ERS events from each other
erd_indx_seppoint = [0, diff(erd_indx)]; % zero-padding for the length
ers_indx_seppoint = [0, diff(ers_indx)]; % zero-padding for the length
erd_sep_indx = find(erd_indx_seppoint~=1);
ers_sep_indx = find(ers_indx_seppoint~=1);

if isempty(erd_sep_indx)
    erd_indx_all = erd_indx;
else
    for i=1:length(erd_sep_indx)
        % these sizes should be variable; DO NOT preallocate
        if i==length(erd_sep_indx)
            erd_indx_all{i} = erd_indx(erd_sep_indx(i):end);
        else
            erd_indx_all{i} = erd_indx(erd_sep_indx(i):erd_sep_indx(i+1));
        end
    end
end
if isempty(ers_sep_indx)
    ers_indx_all = ers_indx;
else
    for i=1:length(ers_sep_indx)
        % this size should be variable; DO NOT preallocate
        if i==length(ers_sep_indx)
            ers_indx_all{i} = ers_indx(ers_sep_indx(i):end);
        else
            ers_indx_all{i} = ers_indx(ers_sep_indx(i):ers_sep_indx(i+1));
        end
    end
end

%% estimating ERD events area enclosed by ERP curve and "mean-confidence_interval"
% values are variable; DO NOT preallocate
for i=1:length(erd_indx_all)
    erd_area_tot{i} = sum((ref_val+intv_tresh(1))*ones(1, length(erd_indx_all{i})));
    erd_area_adj{i} = sum(erp_after(erd_indx_all{i}));
    erd_area{i} = (erd_area_tot{i} - erd_area_adj{i})/(length(erd_indx_all{i})+1); % +1 to prevent NAN
end

%% estimating ERS events area enclosed by ERP curve and "mean+confidence_interval"
% values are variable; DO NOT preallocate
for i=1:length(ers_indx_all)
    ers_area_tot{i} = sum(erp_after(ers_indx_all{i}));
    ers_area_adj{i} = sum((ref_val+intv_tresh(2))*ones(1, length(ers_indx_all{i})));
    ers_area{i} = (ers_area_tot{i} - ers_area_adj{i})/(length(ers_indx_all{i})+1); % +1 to prevent NAN
end

end