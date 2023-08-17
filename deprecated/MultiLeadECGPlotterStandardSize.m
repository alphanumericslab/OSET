function MultiLeadECGPlotterStandardSize(data, ch_names, fs, ref_ch, t1_small, t2_small, Title)
% MultiLeadECGPlotterStandardSize has been deprecated. Use ecg_strip_viewer_standard_grid instead.
warning('MultiLeadECGPlotterStandardSize has been deprecated. Use ecg_strip_viewer_standard_grid instead.');
ecg_strip_viewer_standard_grid(data, ch_names, fs, ref_ch, t1_small, t2_small, Title);
