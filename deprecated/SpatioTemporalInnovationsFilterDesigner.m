function [h, A, S_mean, S_median, S_max, S_min, S_max_min_avg] = SpatioTemporalInnovationsFilterDesigner(template_records, varargin)
% SpatioTemporalInnovationsFilterDesigner has been deprecated. Use spatio_temp_innovs_flt_designer instead.
warning('SpatioTemporalInnovationsFilterDesigner has been deprecated. Use spatio_temp_innovs_flt_designer instead.');
[h, A, S_mean, S_median, S_max, S_min, S_max_min_avg] = spatio_temp_innovs_flt_designer(template_records, varargin{:});
