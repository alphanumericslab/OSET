function [pip, ials, pnn_ss, pnn_as] = hrf(rrn_interval , fs)
%Computes Heart Rate Fragmentation of a NN interval time series.
% HRF for assessing short-term (high-frequency [HF]) heart rate dynamics

% INPUT:
%  rrn_interval: Vector of NN-interval durations millisecond unit
%  fs - (optional input) ECG sampling frequency (default fs=250)
%
% OUTPUT:
% fragmentation metrics:
%
%    PIP - Percentage of inflection points.
%    PSS - Percentage of NN intervals that are in short segments.
%    PAS - Percentage of NN intervals that are in alternation segments of at least 3 intervals.
%    IALS- Inverse average length of segments.
%
%.. [1] Costa, M. D., Davis, R. B., & Goldberger, A. L. (2017). Heart Rate
%   Fragmentation: A New Approach to the Analysis of Cardiac Interbeat Interval
%   Dynamics. Frontiers in Physiology, 8(May), 1–13.
% Author:
%   Sajjad Karimi
%   Emory University, Georgia, USA
%   Email: sajjadkarimi91@gmail.com
%   Date: Mar 18, 2025


if median(rrn_interval) < 5
    rrn_interval = 1000 * rrn_interval;
    warning('input should be in millisecond unit, automaticlly changed from second to millisecond')
end

% r for defining "acceptable" level of noise
if nargin > 1 && ~isempty(fs)
    r = 1000/fs;
else
    r = 4 ; % set to 4 millisecond tolerance by assuming fs = 250
end

%% Computing fragmentation parameters

% Number of RR intervals
N = length(rrn_interval);

% Reshape input into a row vector
rrn_interval = rrn_interval(:)';

% delta NNi: differences of consecutinve RR intervals with adding a dummy
% at first to align length
drri = diff(rrn_interval);
drri = [3*r*sign(drri(1)), drri]; % adding a dummy value

drri(drri<r & drri>-r) = 0; % set small diff values to zero
ind_zero_diff = find(drri==0);
ind_zero_diff(ind_zero_diff==1) =[];

% % nextfill zero ∆NN for correct computing of ∆NNi x ∆NNi-1
% drri(ind_zero_diff) = nan;
% drri = fillmissing(drri,"next"); % if there are two consequtive zeros

% Product of consecutive ∆NN differences to find peaks and valleys
% Inflection points are negative value 
ip_loc = [drri(1:end-1) .* drri(2:end),1] < 0;

% Inflection points are also zero diffs and thier previous neighbors
ip_loc(ind_zero_diff) = true;
ip_loc(ind_zero_diff-1) = true;

% Number of inflection points (where detla NNi changes sign or ∆NN is zero in peaks and valleys)
nip = sum(ip_loc);

% Percentage of inflection points (PIP)
pip = nip / N;


% Indices of inflection points
ip_idx = find(ip_loc);
ind_common = ismember(ind_zero_diff,ip_idx); % index zero ∆NN in set of inflection points

% Length of acceleration/deceleration segments: the difference between inflection point indices
% is the length of the segments. But zero ∆NN should be excluded
segment_lengths = [ip_idx(1), diff(ip_idx)];

segment_lengths(ind_common) = []; % excluding zero ∆NN

% Inverse Average Length of Segments (IALS)
ials = 1 / mean(segment_lengths);

% Number of NN intervals in segments with less than three intervals
short_segment_lengths = segment_lengths(segment_lengths < 3);
nss = sum(short_segment_lengths);

% Percentage of ∆NN intervals in short accelerative and decelerative segments
pnn_ss = nss / N;

diff_pattern = sign(drri); % looking for ABAB or +,-,+,- patterns
ind_alter_candidates = movsum(diff_pattern, [3,0]) == 0 & movsum(abs(diff_pattern), [3,0]) == 4;
index_alters = find(ind_alter_candidates);
ind_alter_candidates(index_alters-1) = true;
ind_alter_candidates(index_alters-2) = true;
ind_alter_candidates(index_alters-3) = true;

% Percentage of NN intervals in alternation segments length > 3 (PAS)
nas = sum(ind_alter_candidates(index_alters-1));
pnn_as = nas / N;

end