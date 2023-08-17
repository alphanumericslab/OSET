function MultiLeadECGPlotter(data, ch_names, fs, ref_ch, t1_small, t2_small, t1_long, t2_long, Title)
% MultiLeadECGPlotter has been deprecated. Use ecg_strip_viewer_multichannel instead.
warning('MultiLeadECGPlotter has been deprecated. Use ecg_strip_viewer_multichannel instead.');
ecg_strip_viewer_multichannel(data, ch_names, fs, ref_ch, t1_small, t2_small, t1_long, t2_long, Title);