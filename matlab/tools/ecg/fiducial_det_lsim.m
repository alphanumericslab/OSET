
function [positions, EXITFLAG] = fiducial_det_lsim(data, ecg_rpeaks_index, fs, varargin)

% ECG fiducial points detector based on LSIM approach
%   [peaks, peak_indexes] = fiducial_det_lsim(data, ecg_rpeaks_index, fs, varargin)
%
%   Inputs:
%       data: Vector of input ECG data
%       ecg_rpeaks_index: R-peak indexes in samples
%       fs: Sampling rate in Hz
%  (optional input):
%        - flag_post_processing: Optional. It is a flag that can be 0 or 1, if it is set to 1 then prioir
%                   information from RR intervals is involved to modifying initail estimated fiducial points by LSIM (default is 1)
%        - flag_prune_P: Optional. It is a flag that can be 0 or 1, if it is set
%                   to 1 then based on a score for p-wave quality assesment the detected p-waves are pruned  (default is 1)
%        - win_qrs: Optional. It is window duration in second unit for extracting diff and std features from ecg signal for QRS onset and offset (default is 0.01 or 10ms)
%        - win_T: Optional. It is window duration in second unit for extracting diff and std features from ecg signal for T-wave onset and offset (default is 0.03 or 30ms)
%        - win_P: Optional. It is window duration in second unit for extracting diff and std features from ecg signal for P-wave onset and offset (default is 0.02 or 20ms)
%        - num_cluster: Number of clusters for FCM (default is 4)
%
%   Outputs:
%	position: struct vector with the detected points locations in samples
%
%   struct description:
%            Pon: P wave onset
%              P: P wave peak
%           Poff: P wave end
%        P_score: An score between [0,1] that indicates the quality of detected P wave
%          QRSon: QRS complex onset
%              R: R wave peak
%         QRSoff: QRS complex end
%            Ton: T wave onset
%              T: first T wave peak
%           Toff: T wave end
% beat_quality_score: an score based on Pearson Correlation indicate the qualty of each ECG beat
%
%  EXITFLAG: It is a struct and the corresponding exit conditions are
%            status: succeeded or failed. In the failed condition, errors can occur in the algorithm when detecting certain positions.
%            message: indicate error message in failed conditions
%   Reference:
%      ........
%
%   Sajjad Karimi, Reza Sameni  2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET
%
%

EXITFLAG.status = 'succeeded';
EXITFLAG.message = [];

if nargin < 3
    error('The first 3 inputs are necessary for the fiducial detection')
end

fs_fd = 1000;
fs_fd = fs;

sample_10ms = round(fs_fd*0.01);
sample_70ms = round(fs_fd*0.07);
sample_100ms = round(fs_fd*0.1);
sample_250ms = round(fs_fd*0.25);
sample_350ms = round(fs_fd*0.35);

% Check optional input arguments
if nargin > 3 && ~isempty(varargin{1})
    flag_post_processing = varargin{1};
else
    flag_post_processing = 0;
end

if nargin > 4 && ~isempty(varargin{2})
    flag_prune_P = varargin{2};
else
    flag_prune_P = 0;
end

if nargin > 5 && ~isempty(varargin{3})
    win_sample_qrs = round(fs*varargin{3});
else
    win_sample_qrs = 2*sample_10ms;
end

if nargin > 6 && ~isempty(varargin{4})
    win_sample_T = round(fs*varargin{4});
else
    win_sample_T = 3*sample_10ms;
end

if nargin > 7 && ~isempty(varargin{5})
    win_sample_P = round(fs*varargin{5});
else
    win_sample_P = 2*sample_10ms;
end

if nargin > 7 && ~isempty(varargin{6})
    num_cluster = varargin{6};
else
    num_cluster = 4;
end


if length(ecg_rpeaks_index)<3
    L = length(ecg_rpeaks_index);

    positions.Pon = nan(1,L);
    positions.P = nan(1,L);
    positions.Poff = nan(1,L);

    positions.QRSon = nan(1,L);
    positions.R = nan(1,L);
    positions.QRSoff = nan(1,L);

    positions.Ton = nan(1,L);
    positions.T = nan(1,L);
    positions.Toff = nan(1,L);

    positions.P_score = zeros(1,L);
    positions.beat_quality_score = zeros(1,L);

    return
end


data = data(:)';
data = data - lp_filter_zero_phase(data, 0.5/fs);
ecg_rpeaks_index = ecg_rpeaks_index(:);
ecg_rpeaks_index_org = ecg_rpeaks_index;
ind_org = 1:length(ecg_rpeaks_index_org);

[beat_quality_score, ecg_rpeaks_index, P_rpeaks_index, max_min_r] = preprocess_rpeaks(data, ecg_rpeaks_index, fs, num_cluster);
ecg_rpeaks_index_org = ecg_rpeaks_index;

if flag_post_processing>0
    index_remove1 = find(beat_quality_score<0.6);
else
    index_remove1 = find(beat_quality_score<0.4);
end
ecg_rpeaks_index(index_remove1) =[];
ind_org(index_remove1) =[];

rr_intervals_ecg = em_interval_calc(ecg_rpeaks_index);
N = min(60,ceil(length(ecg_rpeaks_index)/20));
avg_intervals_ecg = movmean(rr_intervals_ecg,[N,N]);
index_remove = 1+find(((diff(ecg_rpeaks_index(:))./avg_intervals_ecg-1)<-0.4 & ~(diff(ecg_rpeaks_index(:))>0.8*fs)) | diff(ecg_rpeaks_index(:))<0.25*fs);
ecg_rpeaks_index(index_remove) =[];

index_remove2 = ind_org(index_remove);
index_remove = unique([index_remove1(:);index_remove2(:)]);


t_second = (0 : length(data)-1)/fs;


% G = gcd( fs_fd,fs );
% data = resample(data,fs_fd/G,fs/G);

try

    rr_intervals_ecg = em_interval_calc(ecg_rpeaks_index);
    avg_intervals_ecg = movmean(rr_intervals_ecg,[60,60]);
    avg_intervals_ecg = [avg_intervals_ecg;avg_intervals_ecg(end)];
    rr_intervals_ecg = [rr_intervals_ecg;rr_intervals_ecg(end)];


    pqrs_bloks = zeros(length(ecg_rpeaks_index)-2,sample_250ms+sample_70ms+1);
    t_bloks = zeros(length(ecg_rpeaks_index)-2,round(sample_350ms*min(2,max(max(1,avg_intervals_ecg/(2*sample_350ms)))))-sample_70ms);

    for p = 2:length(ecg_rpeaks_index)-1
        this_qrs_index = ecg_rpeaks_index(p)-min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1))):ecg_rpeaks_index(p)+sample_70ms;
        if any(this_qrs_index<1)||any(this_qrs_index>length(data))
            continue;
        end
        pqrs_bloks(p-1,end-length(this_qrs_index)+1:end) = data(this_qrs_index) ;


        this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:ecg_rpeaks_index(p)+min(sample_350ms*min(2,max(1,avg_intervals_ecg(p)/(2*sample_350ms))),floor(0.7*rr_intervals_ecg(p)));
        if any(this_t_index>length(data))
            continue;
        end

        t_bloks(p-1,1:length(this_t_index)) =  data(this_t_index);

    end

    ecg_blocks = [pqrs_bloks,t_bloks];

    max_val = max(prctile(pqrs_bloks(:),99.9),prctile(t_bloks(:),99.9));
    min_val = min(prctile(pqrs_bloks(:),0.1),prctile(t_bloks(:),0.1));

    ecg_blocks(ecg_blocks>max_val) = max_val;
    ecg_blocks(ecg_blocks<min_val) = min_val;
    % ecg_blocks = (ecg_blocks - min_val) / (max_val-min_val);
    for p = 1:size(ecg_blocks,1)
        ecg_blocks(p,:) = ecg_blocks(p,:) - linspace(mean(ecg_blocks(p,1:sample_10ms),2),mean(ecg_blocks(p,end-sample_10ms+1:end),2),size(ecg_blocks,2));
    end

    ecg_blocks(:,1:sample_10ms) = ecg_blocks(:,1:sample_10ms) .* linspace(0.01,1,sample_10ms);
    ecg_blocks(:,end-sample_10ms+1:end) = ecg_blocks(:,end-sample_10ms+1:end) .* linspace(1,0.01,sample_10ms);

    ecg_blocks_normalized = ecg_blocks./sqrt(sum(ecg_blocks.^2,2)); % normalization to unit power
    ecg_blocks_ndiff = [zeros(size(ecg_blocks_normalized,1),win_sample_qrs-ceil(win_sample_qrs/2)), ecg_blocks_normalized(:,win_sample_qrs+1:end) - ecg_blocks_normalized(:,1:end-win_sample_qrs),zeros(size(ecg_blocks_normalized,1),win_sample_qrs-floor(win_sample_qrs/2))];

    fcm_features = ecg_blocks_normalized(:,sample_70ms:end);
    fcm_features = fillmissing(fcm_features,"linear");
    % fcm_features = zscore(fcm_features);

    fcm_options = fcmOptions(NumClusters=num_cluster,MaxNumIteration=25, Exponent=1.1, Verbose=0);
    [fcm_centers, fcm_part_mat] = fcm(fcm_features,fcm_options);
    [~,cluster_fcm] = max(fcm_part_mat);

    num_cls = size(fcm_centers,1);
    index_clustering = [];
    count_cls = zeros(num_cls,1);
    for c = 1:num_cls
        index_temp = find(cluster_fcm==c);
        count_cls(c) = length(index_temp);
        index_clustering = [index_clustering;[index_temp(:),c*ones(length(index_temp),1)]];
    end

    if length(cluster_fcm)<60
        index_clustering(:,2)=1;
    else
        index_clustering_org = index_clustering;
        for c = find(count_cls(:)'<30)
            [~,cls_max] = max(count_cls);
            index_clustering(index_clustering_org(:,2)==c,2) = cls_max;
        end
    end

    num_cls = unique(index_clustering(:,2));

    % check the beat quality Pwave quality
    clear cls_mn_ecg
    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        cls_mn_ecg(c,:) = mean(ecg_blocks(ind_c,sample_70ms:end-sample_100ms),1);
    end


    ecg_denoised_n = 0*data;
    ecg_denoised_ndiff = 0*data;

    pqrs_bloks_n = ecg_blocks_normalized(:,1:size(pqrs_bloks,2));
    t_bloks_n = ecg_blocks_normalized(:,size(pqrs_bloks,2)+1:size(pqrs_bloks,2)+size(t_bloks,2));

    pqrs_bloks_ndiff = ecg_blocks_ndiff(:,1:size(pqrs_bloks,2));
    t_bloks_ndiff = ecg_blocks_ndiff(:,size(pqrs_bloks,2)+1:size(pqrs_bloks,2)+size(t_bloks,2));


    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_rpeaks_index(p)-min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1))):ecg_rpeaks_index(p)+sample_70ms;
        if any(this_qrs_index<1)||any(this_qrs_index>length(data))
            continue;
        end
        ecg_denoised_n(this_qrs_index) = pqrs_bloks_n(p-1,end-length(this_qrs_index)+1:end) ;
        ecg_denoised_ndiff(this_qrs_index) = pqrs_bloks_ndiff(p-1,end-length(this_qrs_index)+1:end) ;

        this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:ecg_rpeaks_index(p)+min(sample_350ms*min(2,max(1,avg_intervals_ecg(p)/(2*sample_350ms))),floor(0.7*rr_intervals_ecg(p)));
        if any(this_t_index>length(data))
            continue;
        end

        ecg_denoised_n(this_t_index) = t_bloks_n(p-1,1:length(this_t_index)) ;
        ecg_denoised_ndiff(this_t_index) = t_bloks_ndiff(p-1,1:length(this_t_index)) ;

    end

    rpeak_intervals = cell(length(ecg_rpeaks_index)-2,1);
    for p = 2:length(ecg_rpeaks_index)-1
        rpeak_intervals{p-1} = (ecg_rpeaks_index(p)-sample_70ms:ecg_rpeaks_index(p)+sample_70ms)';
    end
    rpeak_intervals_mat = cell2mat( rpeak_intervals);
    L = length(rpeak_intervals_mat);

    % data_r = lp_filter_zero_phase(data, 35/fs);
    % data_rh = data_r;
    % data_r = data_r - lp_filter_zero_phase(data_r,0.2/fs);
    % data_hat = [data_r(rpeak_intervals_mat)',ones(L,1)]*(pinv([data_r(rpeak_intervals_mat)',ones(L,1)])*data_rh(rpeak_intervals_mat)');
    % Cxy_base = 1-(norm(data_rh(rpeak_intervals_mat)' - data_hat)/  norm(data_rh(rpeak_intervals_mat)));
    %
    % for f = 1:8
    %     data_r = lp_filter_zero_phase(data, 35/fs);
    %     data_rh = data_r;
    %     data_r = data_r - lp_filter_zero_phase(data_r,f/fs);
    %     data_hat = [data_r(rpeak_intervals_mat)',ones(L,1)]*(pinv([data_r(rpeak_intervals_mat)',ones(L,1)])*data_rh(rpeak_intervals_mat)');
    %
    %     % temp_coef = corrcoef(data_r(:),data(:));
    %     Cxy(f) = 1-(norm(data_rh(rpeak_intervals_mat)' - data_hat)/  norm(data_rh(rpeak_intervals_mat)));
    % end
    % Cxy(f+1) =0;
    % f = find(Cxy<0.9*Cxy_base,1,'first');


    for f = 1:8 % in Hz
        data_r = lp_filter_zero_phase(data, 35/fs);
        data_rh = data_r;
        data_r = data_r - lp_filter_zero_phase(data_r,f/fs);
        data_hat = [data_r(rpeak_intervals_mat)',ones(L,1)]*(pinv([data_r(rpeak_intervals_mat)',ones(L,1)])*data_rh(rpeak_intervals_mat)');

        temp_coef = corrcoef(data_hat,data_rh(rpeak_intervals_mat)');
        % Cxy(f) = 1-(norm(data_rh(rpeak_intervals_mat)' - data_hat)/  norm(data_rh(rpeak_intervals_mat)));
        Cxy(f) = temp_coef(1,2);
    end
    Cxy(f+1) =0;
    f = find(Cxy<0.95,1,'first');
    Cxy(f+1) =0;
    f = find(Cxy<0.95,1,'first');

    data_r = lp_filter_zero_phase(data, 35/fs);
    data_r = data_r - lp_filter_zero_phase(data_r,f/fs);

    ecg_denoised_std = movstd(data_r,[win_sample_qrs,win_sample_qrs]);
    ecg_denoised_ndiff = [zeros(1,win_sample_qrs-ceil(win_sample_qrs/2)), data_r(:,win_sample_qrs+1:end) - data_r(:,1:end-win_sample_qrs),zeros(1,win_sample_qrs-floor(win_sample_qrs/2))];
    ecg_denoised_ndiff = abs(ecg_denoised_ndiff);

    data_r_env = (envelope(data_r,sample_70ms).*ecg_denoised_std/prctile(ecg_denoised_std,90) + ecg_denoised_ndiff)/2;
    data_r_envp = abs(data_r);

    pqrs_bloks_on = cell(2,length(ecg_rpeaks_index)-2);
    pqrs_bloks_on_index = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_on_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_on_2 = cell(1,length(ecg_rpeaks_index)-2);

    first_qrson = nan * ecg_rpeaks_index;
    for p = 2:length(ecg_rpeaks_index)-1

        temp_index = ecg_rpeaks_index(p)-min(2*sample_100ms,round(avg_intervals_ecg(p)/5)):ecg_rpeaks_index(p)-sample_10ms;
        [TF,P] = islocalmax(data_r_envp(temp_index),'MaxNumExtrema',1,'MinProminence',0.5*P_rpeaks_index(p));
        TF = find(TF);

        if isempty(TF)
            gain_qrs = linspace(0.25,1,length(temp_index));
            [TF,P] = islocalmax(data_r_envp(temp_index).*gain_qrs,'MaxNumExtrema',1,'MinProminence',0.25*P_rpeaks_index(p));
            TF = find(TF);
            if ~isempty(TF) && TF<2*length(temp_index)/3
                TF = [];
            end
        end


        if ~isempty(TF) && sign(data_r(temp_index(TF))) ~= sign(data_r(ecg_rpeaks_index(p)))
            starting_peak_p = temp_index(TF);
        else
            starting_peak_p = ecg_rpeaks_index(p);
        end

        this_thr = 0.05*data_r_env(starting_peak_p) + median(data_r_env);
        max_sample_qrson_thr = 2*sample_100ms+1 - find(data_r_env(starting_peak_p-2*sample_100ms:starting_peak_p) < this_thr, 1 ,"last");
        max_sample_qrson_thr = min(max_sample_qrson_thr,sample_100ms+sample_70ms);
        if ~isempty(TF)
            max_sample_qrson_thr = min(max_sample_qrson_thr,sample_70ms);
        else
            max_sample_qrson_thr = min(max_sample_qrson_thr,sample_100ms+sample_70ms);
        end
        if isempty(max_sample_qrson_thr)
            first_qrson(p) = starting_peak_p;
            max_sample_qrson_thr = sample_70ms;
        else
            first_qrson(p) = starting_peak_p-max_sample_qrson_thr;
        end
        max_sample_qrson_thr = 2*max_sample_qrson_thr;

        max_sample_qrson = round((sample_100ms+2*sample_10ms)*max(1,rr_intervals_ecg(p-1)/(2*sample_350ms)));

        max_sample_qrson = max(max_sample_qrson,max_sample_qrson_thr);

        this_qrson_index = starting_peak_p-max_sample_qrson:starting_peak_p;
        if any(this_qrson_index<1)
            pqrs_bloks_on{1,p-1} = ecg_denoised_ndiff(this_qrson_index);
            pqrs_bloks_on{1,p-1} = ecg_denoised_std(this_qrson_index) ;
            pqrs_bloks_on_index{p-1} = this_qrson_index;
            continue;
        end

        temp = ecg_denoised_std(this_qrson_index(1:end-3*sample_10ms));
        if length(temp)>5*sample_10ms
            a = polyfit(1:length(temp),temp,2);
            if a(1)>0
                ind_drop_pq = round(-a(2)/(2*a(1)));
                ind_drop_pq = min(length(temp)-3*sample_10ms, max(1,ind_drop_pq));
            else
                ind_drop_pq =[];
            end
            if isempty(ind_drop_pq)
                ind_drop_pq = 1;
            end

            this_qrson_index = this_qrson_index(ind_drop_pq:end);
        end

        dur_this = round(length(this_qrson_index)/2);
        gain_sig = ones(1,length(this_qrson_index));
        ind_extrm =[];
        rpeak_val = abs(ecg_denoised_n(starting_peak_p));
        if max_min_r(p)<0
            flag_type = 'max';
            [TF,P] = islocalmax(data_r(this_qrson_index(dur_this:end-sample_10ms)),'MaxNumExtrema',1,'MinProminence',0.1*rpeak_val);
            TF = find(TF);
            p_sort = sort(P,"descend");
            if ~isempty(TF) || (p_sort(2)~=0 && p_sort(1)/p_sort(2)>10 && p_sort(1)>0.03*rpeak_val) || (p_sort(2)==0 && p_sort(1)>0.05*rpeak_val)
                if isempty(TF)
                    [~,TF]= max(P);
                end
                ind_extrm = dur_this-1 +TF;
                P = P(TF);
            end
        else
            flag_type = 'min';
            [TF,P] = islocalmin(data_r(this_qrson_index(dur_this:end-sample_10ms)),'MaxNumExtrema',1,'MinProminence',0.1*rpeak_val);
            TF = find(TF);
            p_sort = sort(P,"descend");
            if ~isempty(TF) || (p_sort(2)~=0 && p_sort(1)/p_sort(2)>10 && p_sort(1)>0.03*rpeak_val) || (p_sort(2)==0 && p_sort(1)>0.05*rpeak_val)
                if isempty(TF)
                    [~,TF]= max(P);
                end
                ind_extrm = dur_this-1 +TF;
                P = P(TF);
            end
        end


        if ~isempty(ind_extrm)
            temp = data_r_envp(this_qrson_index(dur_this:end-sample_10ms));
            sample_win_gain = TF - find(temp(1:TF)<temp(TF)/3,1,'last')+1;
            if isempty(sample_win_gain)
                sample_win_gain = sample_10ms;
            end
            temp = min(5,max(0,0.5*(rpeak_val/P)-1))*gausswin(2*sample_win_gain+1)';
            gain_sig(ind_extrm-sample_win_gain:ind_extrm) = gain_sig(ind_extrm-sample_win_gain:ind_extrm)+temp(sample_win_gain+1:end);
        end


        temp = ecg_denoised_std(this_qrson_index);
        thr_pq = mean(temp(end-3*sample_10ms:end))/5;
        [TF,P] = islocalmin(temp(1:end-3*sample_10ms),'MaxNumExtrema',1,'MinProminence',thr_pq/2);
        ind_drop_pq = find(TF);

        if isempty(ind_drop_pq)
            ind_drop_pq = 1;
        else
            ind_drop_pq = find(temp(1:ind_drop_pq)>1.2*temp(ind_drop_pq),1,'last');
            if isempty(ind_drop_pq)
                ind_drop_pq = find(TF);
            end
        end

        this_qrson_index = this_qrson_index(ind_drop_pq:end);
        prior_win = gausswin(4*max(sample_70ms ,length(this_qrson_index))+1)';
        prior_win = prior_win(ceil(length(prior_win)/2)-2*length(this_qrson_index):ceil(length(prior_win)/2)+2*length(this_qrson_index));
        gain_sig = gain_sig(ind_drop_pq:end).*prior_win(length(this_qrson_index)+2:2*length(this_qrson_index)+1);

        pqrs_bloks_on{1,p-1} = data_r_env(this_qrson_index).*gain_sig;
        pqrs_bloks_on{2,p-1} = ecg_denoised_std(this_qrson_index).*gain_sig ;

        pqrs_bloks_on_index{p-1} = this_qrson_index;
        dur_block(p-1,1) = length(this_qrson_index);

        data_r_env(this_qrson_index) = data_r_env(this_qrson_index).*gain_sig;
        ecg_denoised_std(this_qrson_index) = ecg_denoised_std(this_qrson_index).*gain_sig;

        dur_this = round(length(this_qrson_index)/2);
        mean_init_on_1{1,p-1} = [data_r_env(this_qrson_index(1:dur_this));ecg_denoised_std(this_qrson_index(1:dur_this))];
        mean_init_on_2{1,p-1} = [data_r_env(this_qrson_index(dur_this:end));ecg_denoised_std(this_qrson_index(dur_this:end))];

    end

    type_det = 'qrson';
    lock_where = 'last'; % first or last
    obs_seqment_beats = pqrs_bloks_on;
    mean_init_1 = mean_init_on_1;
    mean_init_2 = mean_init_on_2;
    bloks_index = pqrs_bloks_on_index;

    cls_onoff_det = cell(length(num_cls),1);
    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        cls_seqment_beats = obs_seqment_beats(:,ind_c); init_state_1 = mean_init_1(ind_c); init_state_2 = mean_init_2(ind_c);
        v_trans_det = lsim_fit_detect(cls_seqment_beats,type_det, init_state_1, init_state_2);

        if contains(lock_where,'last')
            v_trans_det = dur_block(ind_c) - v_trans_det;
        end

        Md_temp = round(median(v_trans_det(:),'omitmissing'));
        v_trans_det = max(0,50-Md_temp) + v_trans_det;

        v_trans_m = movmedian(mean(v_trans_det,2),[100,100],'omitmissing');
        v_trans_det(:,3) = mean(v_trans_det,2);

        v_trans_filtered = v_trans_det(:,3);
        for p = 1:length(ind_c)
            if abs(v_trans_det(p,1)-v_trans_det(p,2))>sample_10ms
                [val_ms,ind_min] = min(abs(v_trans_det(p,:)-v_trans_m(p)));
                v_trans_filtered(p,1) = v_trans_det(p,ind_min);
            end
        end

        if flag_post_processing>0
            cls_onoff_det{c} = em_interval_calc([zeros(length(v_trans_filtered),1),v_trans_filtered],0.25,0.15,30);
        else
            cls_onoff_det{c} = v_trans_filtered;
        end

        if contains(lock_where,'last')
            cls_onoff_det{c} = dur_block(ind_c) - (cls_onoff_det{c} - max(0,50-Md_temp));
        else
            cls_onoff_det{c} = cls_onoff_det{c} - max(0,50-Md_temp);
        end

        if  contains(type_det , 'on')
            cls_onoff_det{c} = round(cls_onoff_det{c}) ; - ceil(win_sample_qrs/2); % diff correction
        else
            cls_onoff_det{c} = round(cls_onoff_det{c}) ;+ ceil(win_sample_qrs/2); % diff correction
        end
        cls_onoff_det{c}(cls_onoff_det{c}<1) = 1;

    end

    temp = cell2mat(cls_onoff_det);
    temp(index_clustering(:,1)) = temp;
    cls_onoff_det_mat = [nan;temp;nan];

    this_fiducials = nan(1,length(ecg_rpeaks_index));
    for p = 2:length(ecg_rpeaks_index)-1
        this_index = bloks_index{p-1};
        this_fiducials(p) = this_index(min(length(this_index),cls_onoff_det_mat(p)));
        if p==2
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1)>0
                this_fiducials(1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1);
            end
        elseif p==length(ecg_rpeaks_index)-1
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1)<=length(data)
                this_fiducials(p+1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1);
            end
        end
    end

    this_fiducials = min(first_qrson(:)',this_fiducials);
    positions.QRSon = this_fiducials;


    ecg_qrson_index = positions.QRSon(:);
    pqrs_bloks_off = cell(2,length(ecg_rpeaks_index)-2);
    pqrs_bloks_off_index = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_off_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_off_2 = cell(1,length(ecg_rpeaks_index)-2);


    mn_qrson =  mean(ecg_rpeaks_index - ecg_qrson_index, 'omitmissing');
    sample_off_search = round(max(sample_100ms+2*sample_10ms, 3*sample_100ms - mn_qrson));
    first_qrsoff = nan * ecg_rpeaks_index;

    for p = 2:length(ecg_rpeaks_index)-1


        temp_index = ecg_rpeaks_index(p)+sample_10ms:ecg_rpeaks_index(p)+min(2*sample_100ms,round(avg_intervals_ecg(p)/5));

        [TF,P] = islocalmax(data_r_envp(temp_index),'MaxNumExtrema',1,'MinProminence',0.5*P_rpeaks_index(p));
        TF = find(TF);

        if isempty(TF)
            gain_qrs = linspace(1,0.5,length(temp_index));
            [TF,P] = islocalmax(data_r_envp(temp_index).*gain_qrs,'MaxNumExtrema',1,'MinProminence',0.25*P_rpeaks_index(p));
            TF = find(TF);
            P = sort(P,'descend');
            if ~isempty(TF) && TF>length(temp_index)/3
                TF = [];
            end
        end

        change_sign = 1;
        if ~isempty(TF) && sign(data_r(temp_index(TF))) ~= sign(data_r(ecg_rpeaks_index(p)))
            stoping_peak_p = temp_index(TF);
            change_sign = -1;
        else
            stoping_peak_p = ecg_rpeaks_index(p);
        end

        this_thr = 0.05*data_r_env(stoping_peak_p) + median(data_r_env);
        max_sample_qrsoff_thr =  find(data_r_env(stoping_peak_p:stoping_peak_p+2*sample_100ms) < this_thr, 1 ,"first");
        if ~isempty(TF)
            max_sample_qrsoff_thr = min(max_sample_qrsoff_thr,sample_70ms);
        else
            max_sample_qrsoff_thr = min(max_sample_qrsoff_thr,sample_100ms+sample_70ms);
        end
        if isempty(max_sample_qrsoff_thr)
            first_qrsoff(p) = stoping_peak_p;
            max_sample_qrsoff_thr = sample_70ms;
        else
            first_qrsoff(p) = stoping_peak_p + max_sample_qrsoff_thr;
        end

        max_sample_qrsoff_thr = 2*max_sample_qrsoff_thr;
        max_sample_qrsoff = round((sample_100ms+2*sample_10ms)*max(1,rr_intervals_ecg(p-1)/(2*sample_350ms)));
        max_sample_qrsoff = max([max_sample_qrsoff,min(max_sample_qrsoff_thr, sample_off_search)]);

        this_qrsoff_index = stoping_peak_p+sample_10ms:stoping_peak_p + max_sample_qrsoff;
        if any(this_qrsoff_index>length(data))
            continue;
        end

        temp = ecg_denoised_std(this_qrsoff_index);
        if length(temp)>sample_100ms && mean(temp(1:2*sample_10ms)) < 3*mean(temp(end-2*sample_10ms:end))
            a = polyfit(1:length(temp),temp,2);
            if a(1)>0
                ind_drop_st = round(-a(2)/(2*a(1)));
                ind_drop_st = min(length(temp), max(3*sample_10ms,ind_drop_st));
            else
                ind_drop_st =[];
            end
            if isempty(ind_drop_st)
                ind_drop_st = length(temp);
            end

            this_qrsoff_index = this_qrsoff_index(1:ind_drop_st);
        end

        dur_this = round(length(this_qrsoff_index)/2);
        gain_sig = ones(1,length(this_qrsoff_index));
        ind_extrm =[];
        rpeak_val = abs(data(stoping_peak_p));
        if max_min_r(p) <0 && change_sign==1
            flag_type = 'max';

            [TF,P] = islocalmax(data(this_qrsoff_index(sample_10ms:dur_this)),'MaxNumExtrema',1,'MinProminence',0.1*rpeak_val);
            TF = find(TF);
            [p_sort,ind_p] = sort(P,"descend");
            if ~isempty(TF) || (p_sort(2)~=0 && p_sort(1)/p_sort(2)>10 && p_sort(1)>0.03*rpeak_val && ind_p(1)<ind_p(2)) || (p_sort(2)==0 && p_sort(1)>0.05*rpeak_val)
                if isempty(TF)
                    [~,TF]= max(P);
                end
                ind_extrm = sample_10ms-1 +TF;
                P = P(TF);
            end
        elseif change_sign==1

            flag_type = 'min';

            [TF,P] = islocalmin(data(this_qrsoff_index(sample_10ms:dur_this)),'MaxNumExtrema',1,'MinProminence',0.1*rpeak_val);
            TF = find(TF);
            [p_sort,ind_p] = sort(P,"descend");
            if ~isempty(TF) || (p_sort(2)~=0 && p_sort(1)/p_sort(2)>10 && p_sort(1)>0.03*rpeak_val&& ind_p(1)<ind_p(2)) || (p_sort(2)==0 && p_sort(1)>0.05*rpeak_val)
                if isempty(TF)
                    [~,TF]= max(P);
                end
                ind_extrm = sample_10ms-1 +TF;
                P = P(TF);
            end

        end

        if ~isempty(ind_extrm)
            temp = data_r_envp(this_qrsoff_index(sample_10ms:dur_this));
            sample_win_gain = TF - find(temp(1:TF)<temp(TF)/3,1,'first')+1;
            if isempty(sample_win_gain)
                sample_win_gain = sample_10ms;
            end
            temp = min(5,max(0,0.5*(rpeak_val/P)-1))*gausswin(2*sample_win_gain+1)';
            gain_sig(ind_extrm:ind_extrm+sample_win_gain) = gain_sig(ind_extrm:ind_extrm+sample_win_gain)+temp(1:sample_win_gain+1);
        end

        pqrs_bloks_off{1,p-1} = data_r_env(this_qrsoff_index).*gain_sig;
        pqrs_bloks_off{2,p-1} = ecg_denoised_std(this_qrsoff_index).*gain_sig ;
        pqrs_bloks_off_index{p-1} = this_qrsoff_index;
        dur_block(p-1,1) = length(this_qrsoff_index);

        data_r_env(this_qrsoff_index) = data_r_env(this_qrsoff_index).*gain_sig;
        ecg_denoised_std(this_qrsoff_index) = ecg_denoised_std(this_qrsoff_index).*gain_sig;

        dur_this = round(length(this_qrsoff_index)/2);
        mean_init_off_1{1,p-1} = [data_r_env(this_qrsoff_index(1:dur_this));ecg_denoised_std(this_qrsoff_index(1:dur_this))];
        mean_init_off_2{1,p-1} = [data_r_env(this_qrsoff_index(dur_this:end));ecg_denoised_std(this_qrsoff_index(dur_this:end))];

    end

    lock_where = 'first';
    type_det = 'qrsoff';
    obs_seqment_beats = pqrs_bloks_off;
    mean_init_1 = mean_init_off_1;
    mean_init_2 = mean_init_off_2;
    bloks_index = pqrs_bloks_off_index;

    cls_onoff_det = cell(length(num_cls),1);
    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        cls_seqment_beats = obs_seqment_beats(:,ind_c); init_state_1 = mean_init_1(ind_c); init_state_2 = mean_init_2(ind_c);
        v_trans_det = lsim_fit_detect(cls_seqment_beats,type_det, init_state_1, init_state_2);

        if contains(lock_where,'last')
            v_trans_det = dur_block(ind_c) - v_trans_det;
        end

        Md_temp = round(median(v_trans_det(:),'omitmissing'));
        v_trans_det = max(0,50-Md_temp) + v_trans_det;

        v_trans_m = movmedian(mean(v_trans_det,2),[100,100],'omitmissing');
        v_trans_det(:,3) = mean(v_trans_det,2);

        v_trans_filtered = v_trans_det(:,3);
        for p = 1:length(ind_c)
            if abs(v_trans_det(p,1)-v_trans_det(p,2))>sample_10ms
                [val_ms,ind_min] = min(abs(v_trans_det(p,:)-v_trans_m(p)));
                v_trans_filtered(p,1) = v_trans_det(p,ind_min);
            end
        end

        if flag_post_processing>0
            cls_onoff_det{c} = em_interval_calc([zeros(length(v_trans_filtered),1),v_trans_filtered],0.25,0.15,30);
        else
            cls_onoff_det{c} = v_trans_filtered;
        end

        if contains(lock_where,'last')
            cls_onoff_det{c} = dur_block(ind_c) - (cls_onoff_det{c} - max(0,50-Md_temp));
        else
            cls_onoff_det{c} = cls_onoff_det{c} - max(0,50-Md_temp);
        end

        if  contains(type_det , 'on')
            cls_onoff_det{c} = round(cls_onoff_det{c}) ; - ceil(win_sample_qrs/2); % diff correction
        else
            cls_onoff_det{c} = round(cls_onoff_det{c}); + ceil(win_sample_qrs/2); % diff correction
        end
        cls_onoff_det{c}(cls_onoff_det{c}<1) = 1;

    end


    temp = cell2mat(cls_onoff_det);
    temp(index_clustering(:,1)) = temp;
    cls_onoff_det_mat = [nan;temp;nan];

    this_fiducials = nan(1,length(ecg_rpeaks_index));
    for p = 2:length(ecg_rpeaks_index)-1
        this_index = bloks_index{p-1};
        this_fiducials(p) = this_index(min(length(this_index),cls_onoff_det_mat(p)));
        if p==2
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1)>0
                this_fiducials(1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1);
            end
        elseif p==length(ecg_rpeaks_index)-1
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1)<=length(data)
                this_fiducials(p+1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1);
            end
        end
    end

    positions.QRSoff = this_fiducials;

    % ecg_rpeaks_index_p = ecg_rpeaks_index;
    % ecg_rpeaks_index_p(isnan(ecg_rpeaks_index_p)) = [];
    % ecg_QRSoff_index_p = positions.QRSoff ;
    % ecg_QRSoff_index_p(isnan(ecg_QRSoff_index_p)) = [];
    % a = figure('Position', [130 130 1500 800]);
    % lg = {};
    % plot(t_second,data,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
    % hold on
    % plot(t_second(ecg_rpeaks_index_p),data(ecg_rpeaks_index_p),'*',LineWidth=1.5) ;lg = cat(2, lg, {'Rpeaks'});
    % plot(t_second(ecg_QRSoff_index_p),data(ecg_QRSoff_index_p),'o',LineWidth=1.5) ;lg = cat(2, lg, {'Rpeaks'});
    % plot(t_second,data_r,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
    % plot(t_second,ecg_denoised_ndiff,LineWidth=1.5) ;lg = cat(2, lg, {'diff1'});
    % grid on
    % plot(t_second,ecg_denoised_std,LineWidth=1.5) ;lg = cat(2, lg, {'std'});
    % plot(t_second,data_r_env,LineWidth=1.5) ;lg = cat(2, lg, {'std'});

    ecg_qrsoff_index = positions.QRSoff;
    ecg_denoised_nT = data;

    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_qrson_index(p) : ecg_qrsoff_index(p);
        ecg_denoised_nT(this_qrs_index) = linspace(ecg_denoised_nT(this_qrs_index(1)),ecg_denoised_nT(this_qrs_index(end)) , length(this_qrs_index) );

    end

    ecg_denoised_nT = lp_filter_zero_phase(ecg_denoised_nT, 30/fs);

    ecg_T_ndiff = [zeros(1,win_sample_T-ceil(win_sample_T/2)), ecg_denoised_nT(win_sample_T+1:end) - ecg_denoised_nT(1:end-win_sample_T),zeros(1,win_sample_T-floor(win_sample_T/2))];
    ecg_T_ndiff = abs(ecg_T_ndiff);

    std_data = std(ecg_denoised_nT);
    std_diff = mean(ecg_T_ndiff)+2*std(ecg_T_ndiff);
    std_diff2 = mean(ecg_T_ndiff)+1*std(ecg_T_ndiff);
    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_qrson_index(p) : ecg_qrsoff_index(p);
        % abs(ecg_denoised_nT(this_qrs_index(1)))
        if abs(ecg_denoised_nT(this_qrs_index(1))) > std_data/2 && abs(ecg_denoised_nT(this_qrs_index(1)) - ecg_denoised_nT(this_qrs_index(end))) > std_data/2
            [val_max,idx_max] = max(ecg_T_ndiff(ecg_qrson_index(p)-sample_70ms:ecg_qrson_index(p)-1));
            idx_max = sample_70ms - idx_max+1;
            if idx_max<4*sample_10ms && (positions.QRSoff(p) - positions.QRSon(p)) < 3*sample_100ms && (val_max>std_diff || (val_max>std_diff2 && abs(ecg_denoised_nT(this_qrs_index(1))) > std_data && abs(ecg_denoised_nT(this_qrs_index(1)) - ecg_denoised_nT(this_qrs_index(end))) > std_data))
                positions.QRSon(p) = positions.QRSon(p) - 2*idx_max;
            end
        end

    end

    ecg_qrson_index = positions.QRSon;
    ecg_qrsoff_index = positions.QRSoff;
    ecg_denoised_nT = data;

    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_qrson_index(p) : ecg_qrsoff_index(p);
        ecg_denoised_nT(this_qrs_index) = linspace(ecg_denoised_nT(this_qrs_index(1)),ecg_denoised_nT(this_qrs_index(end)) , length(this_qrs_index) );

    end

    ecg_denoised_nT = lp_filter_zero_phase(ecg_denoised_nT, 30/fs);

    ecg_T_ndiff = [zeros(1,win_sample_T-ceil(win_sample_T/2)), ecg_denoised_nT(win_sample_T+1:end) - ecg_denoised_nT(1:end-win_sample_T),zeros(1,win_sample_T-floor(win_sample_T/2))];
    ecg_T_ndiff = abs(ecg_T_ndiff);

    std_data = std(ecg_denoised_nT);
    std_diff = mean(ecg_T_ndiff)+2*std(ecg_T_ndiff);
    std_diff2 = mean(ecg_T_ndiff)+1*std(ecg_T_ndiff);
    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_qrson_index(p) : ecg_qrsoff_index(p);
        if abs(ecg_denoised_nT(this_qrs_index(end))) > std_data/2 && abs(ecg_denoised_nT(this_qrs_index(1)) - ecg_denoised_nT(this_qrs_index(end))) > std_data/2
            [val_max,idx_max] = max(ecg_T_ndiff(ecg_qrsoff_index(p)+1:ecg_qrsoff_index(p)+sample_100ms));
            if idx_max<4*sample_10ms && (positions.QRSoff(p) - positions.QRSon(p)) < 3*sample_100ms && (val_max>std_diff || (val_max>std_diff2 && abs(ecg_denoised_nT(this_qrs_index(end))) > std_data && abs(ecg_denoised_nT(this_qrs_index(1)) - ecg_denoised_nT(this_qrs_index(end))) > std_data))
                positions.QRSoff(p) = positions.QRSoff(p) + 2*idx_max;
            end
        end

    end


    %%  T-wave detection ##################################################################################################################
    %  ================= ##################################################################################################################

    T_bloks_on = cell(2,length(ecg_rpeaks_index)-2);
    T_bloks_on_index = cell(1,length(ecg_rpeaks_index)-2);

    T_bloks_off = cell(2,length(ecg_rpeaks_index)-2);
    T_bloks_off_index = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_Ton_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_Ton_2 = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_Toff_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_Toff_2 = cell(1,length(ecg_rpeaks_index)-2);

    ecg_qrsoff_index = positions.QRSoff;
    ecg_denoised_nT = data;

    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_qrson_index(p) : ecg_qrsoff_index(p);
        ecg_denoised_nT(this_qrs_index) = linspace(ecg_denoised_nT(this_qrs_index(1)),ecg_denoised_nT(this_qrs_index(end)) , length(this_qrs_index) );

    end

    ecg_denoised_nT = lp_filter_zero_phase(ecg_denoised_nT, 10/fs);
    % ecg_denoised_nT = ecg_denoised_nT - lp_filter_zero_phase(ecg_denoised_nT, 0.5/fs);
    % ecg_denoised_nT = sjk_eeg_filter(ecg_denoised_nT, fs_fd,0.5,20);

    ecg_T_std = movstd(ecg_denoised_nT,[win_sample_T,win_sample_T]);
    ecg_T_stds = movstd(ecg_denoised_nT,[2*win_sample_T,2*win_sample_T]);
    ecg_T_ndiff = [zeros(1,win_sample_T-ceil(win_sample_T/2)), ecg_denoised_nT(win_sample_T+1:end) - ecg_denoised_nT(1:end-win_sample_T),zeros(1,win_sample_T-floor(win_sample_T/2))];
    ecg_T_ndiff = abs(ecg_T_ndiff);

    T_peaks = nan(length(ecg_rpeaks_index),1);
    before_peak = nan(length(ecg_rpeaks_index),1);
    after_peak = nan(length(ecg_rpeaks_index),1);
    biphasic_shape = nan(length(ecg_rpeaks_index),1);
    T_type = cell(length(ecg_rpeaks_index),1);

    for p = 2:length(ecg_rpeaks_index)-1

        this_T_index = ecg_qrsoff_index(p)+1:ecg_qrsoff_index(p)+min(sample_350ms*max(1,avg_intervals_ecg(p)/(2*sample_350ms)), ecg_qrson_index(p+1) - ecg_qrsoff_index(p)-sample_70ms);

        if any(this_T_index>length(data))
            continue;
        end

        if length(this_T_index)<2*sample_70ms
            this_T_index = ecg_qrsoff_index(p)+1:ecg_qrsoff_index(p)+sample_250ms;
        end

        temp = ecg_denoised_nT(this_T_index(5*sample_10ms+1:min(length(this_T_index),round(sample_350ms*max(1,avg_intervals_ecg(p)/(3*sample_350ms)) )) ));
        temp2 = ecg_T_stds(this_T_index(5*sample_10ms+1:min(length(this_T_index),round(sample_350ms*max(1,avg_intervals_ecg(p)/(3*sample_350ms)) )) ));

        [TF,Pv_max] = islocalmax(temp,'MaxNumExtrema',2,'MinSeparation',sample_70ms);
        TF_max = find(TF);
        % temp2.*(temp-median(temp(1:20)))/max(temp2)
        [TF,Pv_min] = islocalmin(temp, 'MaxNumExtrema',2,'MinSeparation',sample_70ms);
        TF_min =find(TF);

        if isempty(TF_max) && isempty(TF_min)
            temp(1)=ecg_denoised_nT(this_T_index(1));
            temp = temp  - linspace(temp(1),temp(end),length(temp));

            [TF,Pv_max] = islocalmax(temp,'MaxNumExtrema',2,'MinSeparation',sample_70ms);
            TF_max = find(TF);
            [TF,Pv_min] = islocalmin(temp, 'MaxNumExtrema',2,'MinSeparation',sample_70ms);
            TF_min =find(TF);

        end

        if ~isempty(TF_max)
            P = Pv_max;
            if length(TF_max)==1
                P_max=P(TF_max);
                P(TF_max) = 0;
                [P_maxU, TF_maxU] = max(P);
                TF_max(2) = TF_maxU;
            else
                P_max=P(TF_max(1));
                P_maxU=P(TF_max(2));
                if P_max < P_maxU
                    TF_max = flip(TF_max);
                    P_max=P(TF_max(1));
                    P_maxU=P(TF_max(2));
                end
            end
        else
            P_max = 10^-10;
            P_maxU = 10^-12;
        end

        TF_max = 5*sample_10ms+TF_max;

        if ~isempty(TF_min)
            P = Pv_min;
            if length(TF_min)==1
                P_min=P(TF_min);
                P(TF_min) = 0;
                [P_minU, TF_minU] = max(P);
                TF_min(2) = TF_minU;
            else
                P_min=P(TF_min(1));
                P_minU=P(TF_min(2));
                if P_min < P_minU
                    TF_min = flip(TF_min);
                    P_min=P(TF_min(1));
                    P_minU=P(TF_min(2));
                end
            end

        else
            P_min = 10^-10;
            P_minU = 10^-12;
        end

        TF_min = 5*sample_10ms+TF_min;

        if P_max/P_min > 3
            T_type{p} = 'max';
            before_peak(p) = TF_max(1);
            after_peak(p) = TF_max(1);
            T_peaks(p) = this_T_index(TF_max(1));
        elseif P_min/P_max > 3
            T_type{p} = 'min';
            before_peak(p) = TF_min(1);
            after_peak(p) = TF_min(1);
            T_peaks(p) = this_T_index(TF_min(1));
        elseif (max(P_max/P_min ,P_min/P_max) < 2) && (min(P_min,P_max)/max(P_minU,P_maxU)>10) && (P_max>10^-10 || P_min>10^-10)
            T_type{p} = 'bi-phasic';
            before_peak(p) = min(TF_min(1), TF_max(1));
            after_peak(p) = max(TF_min(1), TF_max(1));
            if P_min/P_max > 1
                biphasic_shape(p) = 1;
                T_peaks(p) = this_T_index(TF_min(1));
            else
                biphasic_shape(p) = 2;
                T_peaks(p) = this_T_index(TF_max(1));
            end
        elseif P_max/P_min > 2 && (TF_max(2)>TF_max(1) && TF_max(2)>TF_min(1) && TF_max(1)<TF_min(1)) % there is U wave
            T_type{p} = 'Umax';
            before_peak(p) = TF_max(1);
            after_peak(p) = TF_max(1);
            T_peaks(p) = this_T_index(TF_max(1));
            biphasic_shape(p) = round(max(min(length(this_T_index),sample_350ms),(TF_max(2)+2*TF_min(1))/3));
        elseif P_min/P_max > 2  && (TF_min(2)>TF_max(1) && TF_min(2)>TF_min(1)&& TF_max(1)>TF_min(1)) % there is U wave
            T_type{p} = 'Umin';
            before_peak(p) = TF_min(1);
            after_peak(p) = TF_min(1);
            T_peaks(p) = this_T_index(TF_min(1));
            biphasic_shape(p) = round(max(min(length(this_T_index),sample_350ms),(TF_min(2)+2*TF_max(1))/3));
        else
            if P_max/P_min >= 1
                T_type{p} = 'Nmax';
                before_peak(p) = TF_max(1);
                after_peak(p) = TF_max(1);
                T_peaks(p) = this_T_index(TF_max(1));
                biphasic_shape(p) = sample_350ms;
            elseif P_min/P_max > 1
                T_type{p} = 'Nmin';
                before_peak(p) = TF_min(1);
                after_peak(p) = TF_min(1);
                T_peaks(p) = this_T_index(TF_min(1));
                biphasic_shape(p) = sample_350ms;
            end
        end

    end

    if flag_post_processing>0
        before_peak_cell = cell(length(num_cls),1);
        after_peak_cell = cell(length(num_cls),1);
        before_peak = before_peak(2:end-1);
        after_peak = after_peak(2:end-1);

        for c = 1:length(num_cls)

            ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
            before_peak_cell{c,1}  = before_peak(ind_c);
            Md_temp = round(median( before_peak_cell{c,1}(:),'omitmissing'));
            before_peak_cell{c,1} = max(0,50-Md_temp) +  before_peak_cell{c,1};
            before_peak_cell{c,1} = em_interval_calc([zeros(length(ind_c),1), before_peak_cell{c,1}],0.25,0.15,60);
            before_peak_cell{c,1} = round( before_peak_cell{c,1})-max(0,50-Md_temp);

            after_peak_cell{c,1}  = after_peak(ind_c);
            Md_temp = round(median(after_peak_cell{c,1}(:),'omitmissing'));
            after_peak_cell{c,1} = max(0,50-Md_temp) + after_peak_cell{c,1};
            after_peak_cell{c,1} = em_interval_calc([zeros(length(ind_c),1),after_peak_cell{c,1}],0.25,0.15,60);
            after_peak_cell{c,1} = round(after_peak_cell{c,1})-max(0,50-Md_temp);

        end

        before_peak = cell2mat(before_peak_cell);
        before_peak(index_clustering(:,1)) = before_peak;
        before_peak = [nan;before_peak;nan];

        after_peak = cell2mat(after_peak_cell);
        after_peak(index_clustering(:,1)) = after_peak;
        after_peak = [nan;after_peak;nan];

    else
        before_peak = round(before_peak);
        after_peak = round(after_peak);
    end

    for p = 2:length(ecg_rpeaks_index)-1

        if isnan(after_peak(p))
            after_peak(p) = 1;
        end
        this_T_index = ecg_qrsoff_index(p)+1:ecg_qrsoff_index(p) + min(ecg_qrson_index(p+1) - ecg_qrsoff_index(p)-sample_70ms,max(sample_350ms*max(1,avg_intervals_ecg(p)/(2*sample_350ms)), after_peak(p)+3*sample_10ms ));

        if any(this_T_index>length(data))
            continue;
        end

        if length(this_T_index)<2*sample_70ms
            this_T_index = ecg_qrsoff_index(p)+1:ecg_qrsoff_index(p)+sample_250ms;
        end

        T_dur = length(this_T_index);
        if contains(T_type{p},'max')
            T_peaks(p) = this_T_index(before_peak(p));
        elseif contains(T_type{p},'min')
            T_peaks(p) = this_T_index(before_peak(p));
        elseif contains(T_type{p},'bi-phasic')
            T_peaks(p) = this_T_index(before_peak(p));
        elseif contains(T_type{p},'Umax') ||  contains(T_type{p},'Umin')
            T_peaks(p) = this_T_index(before_peak(p));
            this_T_index = this_T_index(1: min(T_dur,max(after_peak(p)+sample_70ms,biphasic_shape(p))));
        elseif contains(T_type{p},'Nmax') ||  contains(T_type{p},'Nmin')
            T_peaks(p) = this_T_index(before_peak(p));
        end

        try
            this_bT_index = this_T_index(1:before_peak(p)-sample_10ms);
        catch
            this_bT_index = this_T_index(1:before_peak(p)/2);
        end

        dur_this = ceil(length(this_bT_index)/3);
        thr_ST = max(ecg_T_std(this_bT_index(2*dur_this:end)))/2;
        ind_drop_st = min(dur_this,find(ecg_T_std(this_bT_index)<thr_ST,1,"first"));
        if isempty(ind_drop_st)
            ind_drop_st = sample_10ms;
        end

        this_bT_index = this_bT_index(ind_drop_st:end);


        T_bloks_on{1,p-1} = ecg_T_ndiff(this_bT_index);
        T_bloks_on{2,p-1} = ecg_T_std(this_bT_index);

        T_bloks_on_index{p-1} = this_bT_index;

        dur_this = ceil(length(this_bT_index)/3);

        mean_init_Ton_1{1,p-1} = [ecg_T_ndiff(this_bT_index(1:dur_this));ecg_T_std(this_bT_index(1:dur_this))];
        mean_init_Ton_2{1,p-1} = [ecg_T_ndiff(this_bT_index(dur_this:end));ecg_T_std(this_bT_index(dur_this:end))];

        this_aT_index = this_T_index(after_peak(p)+win_sample_T:end);

        if isempty(this_aT_index)
            this_aT_index = this_T_index(ceil(before_peak(p)/2):end);
        elseif length(this_aT_index)<2*sample_10ms
            this_aT_index = this_T_index(after_peak(p):end);
        end


        temp = ecg_T_std(this_aT_index);
        if length(temp)>5*sample_10ms
            thr_tp = mean(temp(1:3*sample_10ms))/5;
            [TF,P] = islocalmin(temp(3*sample_10ms+1:end),'MaxNumExtrema',1,'MinProminence',thr_tp/2);
            ind_drop_pq = 3*sample_10ms+find(TF);

            if isempty(ind_drop_pq)
                ind_drop_pq = length(temp);
            else
                ind_drop_pq = ind_drop_pq + find(temp(ind_drop_pq+1:end)>1.2*temp(ind_drop_pq),1,'first');
                if isempty(ind_drop_pq)
                    ind_drop_pq = 3*sample_10ms+find(TF);
                end
            end

            this_aT_index = this_aT_index(1:ind_drop_pq);
        end

        prior_win = gausswin(4*max(sample_70ms ,length(this_aT_index))+1)';
        prior_win = prior_win(ceil(length(prior_win)/2)-2*length(this_aT_index):ceil(length(prior_win)/2)+2*length(this_aT_index));
        gain_sig =prior_win(2*length(this_aT_index)+1:3*length(this_aT_index));

        T_bloks_off{1,p-1} = ecg_T_ndiff(this_aT_index).*gain_sig;
        T_bloks_off{2,p-1} = ecg_T_std(this_aT_index).*gain_sig;

        T_bloks_off_index{p-1} = this_aT_index;

        ecg_T_ndiff(this_aT_index) = ecg_T_ndiff(this_aT_index).*gain_sig;
        ecg_T_std(this_aT_index) = ecg_T_std(this_aT_index).*gain_sig;

        dur_this = ceil(length(this_aT_index)/2);
        mean_init_Toff_1{1,p-1} = [ecg_T_ndiff(this_aT_index(1:dur_this));ecg_T_std(this_aT_index(1:dur_this))];
        mean_init_Toff_2{1,p-1} = [ecg_T_ndiff(this_aT_index(dur_this:end));ecg_T_std(this_aT_index(dur_this:end))];


    end


    if T_peaks(2)- ecg_rpeaks_index(2) + ecg_rpeaks_index(1)>0
        T_peaks(1) = T_peaks(2)- ecg_rpeaks_index(2) + ecg_rpeaks_index(1);
    end

    if T_peaks(end-1)- ecg_rpeaks_index(end-1) + ecg_rpeaks_index(end)<=length(data)
        T_peaks(end) = T_peaks(end-1)- ecg_rpeaks_index(end-1) + ecg_rpeaks_index(end);
    end

    positions.T = T_peaks';

    lock_where = 'first';
    type_det = 'on';
    obs_seqment_beats = T_bloks_on;
    mean_init_1 = mean_init_Ton_1;
    mean_init_2 = mean_init_Ton_2;
    bloks_index = T_bloks_on_index;

    cls_onoff_det = cell(length(num_cls),1);
    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        cls_seqment_beats = obs_seqment_beats(:,ind_c); init_state_1 = mean_init_1(ind_c); init_state_2 = mean_init_2(ind_c);
        v_trans_det = lsim_fit_detect(cls_seqment_beats,type_det, init_state_1, init_state_2);

        if contains(lock_where,'last')
            v_trans_det = dur_block(ind_c) - v_trans_det;
        end

        Md_temp = round(median(v_trans_det(:),'omitmissing'));
        v_trans_det = max(0,50-Md_temp) + v_trans_det;

        v_trans_m = movmedian(mean(v_trans_det,2),[100,100],'omitmissing');
        v_trans_det(:,3) = mean(v_trans_det,2);

        v_trans_filtered = v_trans_det(:,3);
        for p = 1:length(ind_c)
            if abs(v_trans_det(p,1)-v_trans_det(p,2))>sample_10ms
                [val_ms,ind_min] = min(abs(v_trans_det(p,:)-v_trans_m(p)));
                v_trans_filtered(p,1) = v_trans_det(p,ind_min);
            end
        end

        if flag_post_processing>0
            cls_onoff_det{c} = em_interval_calc([zeros(length(v_trans_filtered),1),v_trans_filtered],0.25,0.15,30);
        else
            cls_onoff_det{c} = v_trans_filtered;
        end

        if contains(lock_where,'last')
            cls_onoff_det{c} = dur_block(ind_c) - (cls_onoff_det{c} - max(0,50-Md_temp));
        else
            cls_onoff_det{c} = cls_onoff_det{c} - max(0,50-Md_temp);
        end

        if  contains(type_det , 'on')
            cls_onoff_det{c} = round(cls_onoff_det{c}) ;- ceil(win_sample_T/2); % diff correction
        else
            cls_onoff_det{c} = round(cls_onoff_det{c}) ;+ ceil(win_sample_T/2); % diff correction
        end

        cls_onoff_det{c}(cls_onoff_det{c}<1) = 1;

    end


    temp = cell2mat(cls_onoff_det);
    temp(index_clustering(:,1)) = temp;
    cls_onoff_det_mat = [nan;temp;nan];

    this_fiducials = nan(1,length(ecg_rpeaks_index));
    for p = 2:length(ecg_rpeaks_index)-1
        this_index = bloks_index{p-1};
        this_fiducials(p) = this_index(min(length(this_index),cls_onoff_det_mat(p)));
        if p==2
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1)>0
                this_fiducials(1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1);
            end
        elseif p==length(ecg_rpeaks_index)-1
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1)<=length(data)
                this_fiducials(p+1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1);
            end
        end
    end

    positions.Ton = this_fiducials;

    lock_where = 'first';
    type_det = 'off';
    obs_seqment_beats = T_bloks_off;
    mean_init_1 = mean_init_Toff_1;
    mean_init_2 = mean_init_Toff_2;
    bloks_index = T_bloks_off_index;

    cls_onoff_det = cell(length(num_cls),1);
    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        cls_seqment_beats = obs_seqment_beats(:,ind_c); init_state_1 = mean_init_1(ind_c); init_state_2 = mean_init_2(ind_c);
        v_trans_det = lsim_fit_detect(cls_seqment_beats,type_det, init_state_1, init_state_2);

        if contains(lock_where,'last')
            v_trans_det = dur_block(ind_c) - v_trans_det;
        end

        Md_temp = round(median(v_trans_det(:),'omitmissing'));
        v_trans_det = max(0,50-Md_temp) + v_trans_det;

        v_trans_m = movmedian(mean(v_trans_det,2),[100,100],'omitmissing');
        v_trans_det(:,3) = mean(v_trans_det,2);

        v_trans_filtered = v_trans_det(:,3);
        for p = 1:length(ind_c)
            if abs(v_trans_det(p,1)-v_trans_det(p,2))>sample_10ms
                [val_ms,ind_min] = min(abs(v_trans_det(p,:)-v_trans_m(p)));
                v_trans_filtered(p,1) = v_trans_det(p,ind_min);
            end
        end

        if flag_post_processing>0
            cls_onoff_det{c} = em_interval_calc([zeros(length(v_trans_filtered),1),v_trans_filtered],0.25,0.15,30);
        else
            cls_onoff_det{c} = v_trans_filtered;
        end

        if contains(lock_where,'last')
            cls_onoff_det{c} = dur_block(ind_c) - (cls_onoff_det{c} - max(0,50-Md_temp));
        else
            cls_onoff_det{c} = cls_onoff_det{c} - max(0,50-Md_temp);
        end

        if  contains(type_det , 'on')
            cls_onoff_det{c} = round(cls_onoff_det{c}) ;- ceil(win_sample_T/2); % diff correction
        else
            cls_onoff_det{c} = round(cls_onoff_det{c}) ;  % diff correction
        end

        cls_onoff_det{c}(cls_onoff_det{c}<1) = 1;

    end


    temp = cell2mat(cls_onoff_det);
    temp(index_clustering(:,1)) = temp;
    cls_onoff_det_mat = [nan;temp;nan];

    this_fiducials = nan(1,length(ecg_rpeaks_index));
    for p = 2:length(ecg_rpeaks_index)-1
        this_index = bloks_index{p-1};
        this_fiducials(p) = this_index(min(length(this_index),cls_onoff_det_mat(p)));
        if p==2
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1)>0
                this_fiducials(1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1);
            end
        elseif p==length(ecg_rpeaks_index)-1
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1)<=length(data)
                this_fiducials(p+1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1);
            end
        end
    end

    positions.Toff = this_fiducials;


    % ecg_QRSoff_index_p = positions.QRSoff ;
    % ecg_QRSoff_index_p(isnan(ecg_QRSoff_index_p)) = [];
    % a = figure('Position', [130 130 1500 800]);
    % lg = {};
    % plot(t_second,data,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
    % hold on
    % plot(t_second,ecg_denoised_nT,LineWidth=1.5) ;lg = cat(2, lg, {'ECG-T'});
    % plot(t_second(ecg_QRSoff_index_p),ecg_denoised_nT(ecg_QRSoff_index_p),'*',LineWidth=1.5) ;lg = cat(2, lg, {'ECG-T'});
    % plot(t_second,ecg_T_ndiff,LineWidth=1.5) ;lg = cat(2, lg, {'diff1'});
    % grid on
    % plot(t_second,ecg_T_std,LineWidth=1.5) ;lg = cat(2, lg, {'std'});

    %% P-wave Detection

    ecg_Toff_index = positions.Toff;

    ecg_denoised_nP = nan*data;
    for p = 2:length(ecg_rpeaks_index)-1

        temp_back = min( 2*sample_250ms, max(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
        temp_back = min(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms , temp_back);
        if temp_back<=sample_70ms
            temp_back = sample_100ms;
        end

        this_P_index = ecg_qrson_index(p)-temp_back: ecg_qrson_index(p)-sample_10ms;
        if any(this_P_index<1)||any(this_P_index>length(data))
            continue;
        end

        ecg_denoised_nP(this_P_index) = data(this_P_index) ;
    end

    base_value = median(ecg_denoised_nP,'omitmissing');

    ecg_denoised_nP = fillmissing(ecg_denoised_nP,'linear');
    ecg_denoised_nP = lp_filter_zero_phase(ecg_denoised_nP, 20/fs);
    % ecg_denoised_nP = ecg_denoised_nP - lp_filter_zero_phase(ecg_denoised_nP, 0.1/fs);

    ecg_P_std = movstd(ecg_denoised_nP,[win_sample_P,win_sample_P]);
    ecg_P_ndiff = [zeros(1,win_sample_P-ceil(win_sample_P/2)), ecg_denoised_nP(win_sample_P+1:end) - ecg_denoised_nP(1:end-win_sample_P),zeros(1,win_sample_P-floor(win_sample_P/2))];
    ecg_P_ndiff = abs(ecg_P_ndiff);

    ecg_P_env = (envelope(ecg_denoised_nP,sample_70ms).*ecg_P_std/prctile(ecg_P_std,90) + ecg_P_ndiff)/2;

    % a = figure('Position', [130 130 1500 800]);
    % lg = {};
    % % plot(t_second,data/3,LineWidth=1) ;lg = cat(2, lg, {'ECG'});
    % plot(t_second,ecg_denoised_nP,LineWidth=1.5) ;lg = cat(2, lg, {'ECG-P'});
    % hold on
    % plot(t_second,ecg_P_ndiff,LineWidth=1.5) ;lg = cat(2, lg, {'diff1'});
    % grid on
    % plot(t_second,ecg_P_std,LineWidth=1.5) ;lg = cat(2, lg, {'std'});
    % plot(t_second,ecg_P_env,LineWidth=1.5) ;lg = cat(2, lg, {'env+diff'});

    P_bloks_off = cell(2,length(ecg_rpeaks_index)-2);
    P_bloks_off_index = cell(1,length(ecg_rpeaks_index)-2);

    P_bloks_on = cell(2,length(ecg_rpeaks_index)-2);
    P_bloks_on_index = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_Poff_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_Poff_2 = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_Pon_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_Pon_2 = cell(1,length(ecg_rpeaks_index)-2);


    % ecg_denoised_nP = 0*data;

    max_indexes = nan(length(ecg_rpeaks_index)-2,1);
    max_P_vals = nan(length(ecg_rpeaks_index)-2,1);
    min_indexes = nan(length(ecg_rpeaks_index)-2,1);
    min_P_vals = nan(length(ecg_rpeaks_index)-2,1);
    P_wave_blocks = zeros(length(ecg_rpeaks_index)-2, 2*sample_250ms);
    P_wave_blocks_std = zeros(length(ecg_rpeaks_index)-2, 2*sample_250ms);

    for p = 2:length(ecg_rpeaks_index)-1

        temp_back = min( 2*sample_250ms, max(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
        temp_back = min(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms , temp_back);
        if temp_back<=sample_70ms
            temp_back = sample_100ms;
        end

        this_P_index = ecg_qrson_index(p)-temp_back: ecg_qrson_index(p)-sample_10ms;
        if any(this_P_index<1)||any(this_P_index>length(data))
            continue;
        end

        % ecg_denoised_nP(this_P_index) = data(this_P_index) - linspace(data(this_P_index(1)),data(this_P_index(end)) , length(this_P_index) );
        temp_pindex = ecg_qrson_index(p)-2*sample_250ms-sample_10ms+1:ecg_qrson_index(p)-sample_10ms;
        temp_pindex(temp_pindex<1)=[];
        % temp_pindex(temp_pindex>length(data))=[];

        P_wave_blocks(p-1,end-length(temp_pindex)+1:end) = ecg_P_env(temp_pindex);
        P_wave_blocks_std(p-1,end-length(temp_pindex)+1:end) = ecg_P_std(temp_pindex);
        % ecg_denoised_nP(this_P_index) = ecg_denoised_nP(this_P_index)/norm(ecg_denoised_nP(this_P_index));

        this_P_index = flip(this_P_index);

        [TF,P] = islocalmax(ecg_denoised_nP(this_P_index),'MaxNumExtrema',1);
        TF_max = find(TF);
        if ~isempty(TF_max)
            max_indexes(p-1) = TF_max;
            max_P_vals(p-1) = P(TF_max);
        end

        [TF,P] = islocalmin(ecg_denoised_nP(this_P_index),'MaxNumExtrema',1);
        TF_min = find(TF);
        if ~isempty(TF_min)
            min_indexes(p-1) = TF_min;
            min_P_vals(p-1) = P(TF_min);
        end

    end



    noisy_pwave = 0;
    clear mn_std_index

    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);

        if abs(mean(max_P_vals(ind_c),'omitmissing')) > abs(mean(min_P_vals(ind_c),'omitmissing'))
            std_index = max_indexes;
        else
            std_index = min_indexes;
        end

        mn_std_index(c) = mean(movstd(std_index,[60,60]),'omitmissing');

        if mean(mn_std_index(c),'omitmissing')>sample_10ms
            noisy_pwave(c) = 1;
            P_wave_blocks(ind_c,:) = movmean(P_wave_blocks(ind_c,:),[30,30]);
            P_wave_blocks_std(ind_c,:) = movmean(P_wave_blocks_std(ind_c,:),[30,30]);
        else
            noisy_pwave(c) = 0;
        end

    end


    for p = 2:length(ecg_rpeaks_index)-1

        temp_back = min( 2*sample_250ms, max(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
        temp_back = min(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms , temp_back);
        if temp_back<=sample_70ms
            temp_back = sample_100ms;
        end
        this_P_index = ecg_qrson_index(p)-temp_back: ecg_qrson_index(p)-sample_10ms;
        if any(this_P_index<1)||any(this_P_index>length(data))
            continue;
        end

        ecg_P_env(this_P_index) = P_wave_blocks(p-1,end-length(this_P_index)+1:end) ;
        ecg_P_std(this_P_index) = P_wave_blocks_std(p-1,end-length(this_P_index)+1:end);
    end




    P_peaks = nan(length(ecg_rpeaks_index),1);
    Pbefore_peak = nan(length(ecg_rpeaks_index),1);
    Pafter_peak = nan(length(ecg_rpeaks_index),1);
    Pbiphasic_shape = nan(length(ecg_rpeaks_index),1);
    P_type = cell(length(ecg_rpeaks_index),1);
    P_type{1} = 'x';
    P_type{end} = 'x';

    for p = 2:length(ecg_rpeaks_index)-1

        % this_P_index = ecg_qrsoff_index(p)-round(min(0.66*(ecg_Toff_index(p-1)-ecg_qrsoff_index(p))),  (sample_100ms+sample_70ms)*max(1,avg_intervals_ecg(p)/(2*sample_350ms))):ecg_qrsoff_index(p);

        temp_back = min( 2*sample_250ms, max(min(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms,floor(0.4*rr_intervals_ecg(p-1))) ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
        temp_back = min(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms , temp_back);
        if temp_back<=sample_70ms
            temp_back = sample_100ms;
        end
        this_P_index = ecg_qrson_index(p)-temp_back: ecg_qrson_index(p);

        if any(this_P_index>length(data))
            continue;
        end

        if length(this_P_index)<2*sample_10ms
            this_P_index = ecg_qrson_index(p)-sample_100ms-sample_70ms: ecg_qrson_index(p);
        end

        gain_P = linspace(1,0.5,length(this_P_index));
        this_P_index = flip(this_P_index);
        [TF,P] = islocalmax(ecg_denoised_nP(this_P_index),'MaxNumExtrema',2,'MinSeparation',sample_10ms);
        P = P.*gain_P;
        TF_max = find(TF);
        if ~isempty(TF_max)
            if length(TF_max)==1
                P_max=P(TF_max);
                P(TF_max) = 0;
                [P_maxU, TF_maxU] = max(P);
                TF_max(2) = TF_maxU;
            else
                P_max=P(TF_max(1));
                P_maxU=P(TF_max(2));
                if P_max < P_maxU
                    TF_max = flip(TF_max);
                    P_max=P(TF_max(1));
                    P_maxU=P(TF_max(2));
                end
            end
        else
            P_max = 10^-10;
            P_maxU = 10^-12;
            TF_max = ceil(length(this_P_index)/2);
        end

        [TF,P] = islocalmin(ecg_denoised_nP(this_P_index),'MaxNumExtrema',2,'MinSeparation',sample_10ms);
        P = P.*gain_P;
        TF_min = find(TF);
        if ~isempty(TF_min)
            if length(TF_min)==1
                P_min=P(TF_min);
                P(TF_min) = 0;
                [P_minU, TF_minU] = max(P);
                TF_min(2) = TF_minU;
            else
                P_min=P(TF_min(1));
                P_minU=P(TF_min(2));
                if P_min < P_minU
                    TF_min = flip(TF_min);
                    P_min=P(TF_min(1));
                    P_minU=P(TF_min(2));
                end
            end

        else
            P_min = 10^-10;
            P_minU = 10^-12;
            TF_min = ceil(length(this_P_index)/2);
        end

        if P_max/P_min >= 1
            P_type{p} = 'max';
            Pbefore_peak(p) = TF_max(1);
            Pafter_peak(p) = TF_max(1);
            P_peaks(p) = this_P_index(TF_max(1));
        elseif P_min/P_max > 1
            P_type{p} = 'min';
            Pbefore_peak(p) = TF_min(1);
            Pafter_peak(p) = TF_min(1);
            P_peaks(p) = this_P_index(TF_min(1));
        end

    end

    if flag_post_processing>0
        before_peak_cell = cell(length(num_cls),1);
        after_peak_cell = cell(length(num_cls),1);
        Pbefore_peak = Pbefore_peak(2:end-1);
        Pafter_peak = Pafter_peak(2:end-1);

        for c = 1:length(num_cls)

            ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
            before_peak_cell{c,1}  = Pbefore_peak(ind_c);
            Md_temp = round(median( before_peak_cell{c,1}(:),'omitmissing'));
            before_peak_cell{c,1} = max(0,50-Md_temp) +  before_peak_cell{c,1};
            before_peak_cell{c,1} = em_interval_calc([zeros(length(ind_c),1), before_peak_cell{c,1}],0.15,0.10,60);
            before_peak_cell{c,1} = round( before_peak_cell{c,1})-max(0,50-Md_temp);

            after_peak_cell{c,1}  = Pafter_peak(ind_c);
            Md_temp = round(median(after_peak_cell{c,1}(:),'omitmissing'));
            after_peak_cell{c,1} = max(0,50-Md_temp) + after_peak_cell{c,1};
            after_peak_cell{c,1} = em_interval_calc([zeros(length(ind_c),1),after_peak_cell{c,1}],0.15,0.10,60);
            after_peak_cell{c,1} = round(after_peak_cell{c,1})-max(0,50-Md_temp);

        end

        Pbefore_peak = cell2mat(before_peak_cell);
        Pbefore_peak(index_clustering(:,1)) = Pbefore_peak;
        Pbefore_peak = [nan;Pbefore_peak;nan];

        Pafter_peak = cell2mat(after_peak_cell);
        Pafter_peak(index_clustering(:,1)) = Pafter_peak;
        Pafter_peak = [nan;Pafter_peak;nan];

    end



    for p = 2:length(ecg_rpeaks_index)-1

        if isnan(Pafter_peak(p))
            Pafter_peak(p) = 1;
        end

        temp_back = min( 2*sample_250ms, max(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
        temp_back = min(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms , temp_back);
        if temp_back<=sample_70ms
            temp_back = sample_100ms;
        end

        this_P_index = ecg_qrson_index(p)-temp_back: ecg_qrson_index(p);

        if any(this_P_index>length(data))
            continue;
        end

        if length(this_P_index)<2*sample_10ms
            this_P_index = ecg_qrson_index(p)-sample_100ms-sample_70ms: ecg_qrson_index(p);
        end

        this_P_index = flip(this_P_index);

        T_dur = length(this_P_index);
        if contains(P_type{p},'max')
            P_peaks(p) = this_P_index(min(T_dur,Pbefore_peak(p)));
        elseif contains(P_type{p},'min')
            P_peaks(p) = this_P_index(min(T_dur,Pbefore_peak(p)));
        end


        try
            this_bP_index = this_P_index(1:Pbefore_peak(p));
        catch
            this_bP_index = this_P_index(1:round(Pbefore_peak(p)/2));
        end

        dur_this = ceil(length(this_bP_index)/3);
        thr_QP = max(ecg_P_std(this_bP_index(2*dur_this:end)))/5;
        ind_drop_qp = min(dur_this,find(ecg_P_std(this_bP_index)<thr_QP,1,"first"));
        if isempty(ind_drop_qp)
            ind_drop_qp = 1;
        end

        this_bP_index = this_bP_index(ind_drop_qp:end);

        P_bloks_off{1,p-1} = ecg_P_env(this_bP_index); %P_bloks_on{1,p-1}(:,1:2*sample_10ms) = P_bloks_on{1,p-1}(:,1:2*sample_10ms).*linspace(0,1,2*sample_10ms);
        P_bloks_off{2,p-1} = ecg_P_std(this_bP_index);  %P_bloks_on{2,p-1}(:,1:2*sample_10ms) = P_bloks_on{2,p-1}(:,1:2*sample_10ms).*linspace(0,1,2*sample_10ms);
        P_bloks_off_index{p-1} = this_bP_index;

        dur_this = ceil(length(this_bP_index)/3);
        mean_init_Poff_1{1,p-1} = [ecg_P_env(this_bP_index(1:dur_this));ecg_P_std(this_bP_index(1:dur_this))];
        mean_init_Poff_2{1,p-1} = [ecg_P_env(this_bP_index(dur_this:end));ecg_P_std(this_bP_index(dur_this:end))];

        this_aP_index = this_P_index(Pafter_peak(p):end);

        if isempty(this_aP_index)
            this_aP_index = this_P_index(ceil(Pbefore_peak(p)/2):end);
        elseif length(this_aP_index)<2*sample_10ms
            this_aP_index = this_P_index(Pafter_peak(p)-sample_10ms:end);
        end

        temp = ecg_P_std(this_aP_index);
        if length(temp)>5*sample_10ms
            thr_tp = mean(temp(1:2*sample_10ms))/10;
            [TF,P] = islocalmin(temp(2*sample_10ms+1:end),'MaxNumExtrema',1,'MinProminence',thr_tp);
            ind_drop_pq = 2*sample_10ms+find(TF);

            if isempty(ind_drop_pq)
                ind_drop_pq = length(temp);
            else
                ind_drop_pq = ind_drop_pq + find(temp(ind_drop_pq+1:end)>1.2*temp(ind_drop_pq),1,'first');
                if isempty(ind_drop_pq)
                    ind_drop_pq = 2*sample_10ms+find(TF);
                end
            end

            this_aP_index = this_aP_index(1:ind_drop_pq);
        end

        prior_win = gausswin(4*max(sample_70ms ,length(this_aP_index))+1)';
        prior_win = prior_win(ceil(length(prior_win)/2)-2*length(this_aP_index):ceil(length(prior_win)/2)+2*length(this_aP_index));
        gain_sig =prior_win(2*length(this_aP_index)+1:3*length(this_aP_index));

        P_bloks_on{1,p-1} = ecg_P_env(this_aP_index).*gain_sig;
        P_bloks_on{2,p-1} = ecg_P_std(this_aP_index).*gain_sig;
        P_bloks_on_index{p-1} = this_aP_index;


        ecg_P_env(this_aP_index) = ecg_P_env(this_aP_index).*gain_sig;
        ecg_P_std(this_aP_index) = ecg_P_std(this_aP_index).*gain_sig;

        dur_this = ceil(length(this_aP_index)/2);
        mean_init_Pon_1{1,p-1} = [ecg_P_env(this_aP_index(1:dur_this));ecg_P_std(this_aP_index(1:dur_this))];
        mean_init_Pon_2{1,p-1} = [ecg_P_env(this_aP_index(dur_this:end));ecg_P_std(this_aP_index(dur_this:end))];


    end



    if P_peaks(2)- ecg_rpeaks_index(2) + ecg_rpeaks_index(1)>0
        P_peaks(1) = P_peaks(2)- ecg_rpeaks_index(2) + ecg_rpeaks_index(1);
    end

    if P_peaks(end-1)- ecg_rpeaks_index(end-1) + ecg_rpeaks_index(end)<=length(data)
        P_peaks(end) = P_peaks(end-1)- ecg_rpeaks_index(end-1) + ecg_rpeaks_index(end);
    end

    positions.P = P_peaks';

    lock_where = 'first';
    type_det = 'on';
    obs_seqment_beats = P_bloks_off;
    mean_init_1 = mean_init_Poff_1;
    mean_init_2 = mean_init_Poff_2;
    bloks_index = P_bloks_off_index;

    cls_onoff_det = cell(length(num_cls),1);
    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        cls_seqment_beats = obs_seqment_beats(:,ind_c); init_state_1 = mean_init_1(ind_c); init_state_2 = mean_init_2(ind_c);
        v_trans_det = lsim_fit_detect(cls_seqment_beats,type_det, init_state_1, init_state_2);

        if noisy_pwave(c)>0
            v_trans_det(:,2) = v_trans_det(:,1);
        end

        if contains(lock_where,'last')
            v_trans_det = dur_block(ind_c) - v_trans_det;
        end

        Md_temp = round(median(v_trans_det(:),'omitmissing'));
        v_trans_det = max(0,50-Md_temp) + v_trans_det;

        v_trans_m = movmedian(mean(v_trans_det,2),[100,100],'omitmissing');
        v_trans_det(:,3) = mean(v_trans_det,2);

        v_trans_filtered = v_trans_det(:,3);
        for p = 1:length(ind_c)
            if abs(v_trans_det(p,1)-v_trans_det(p,2))>sample_10ms
                [val_ms,ind_min] = min(abs(v_trans_det(p,:)-v_trans_m(p)));
                v_trans_filtered(p,1) = v_trans_det(p,ind_min);
            end
        end

        if flag_post_processing>0
            cls_onoff_det{c} = em_interval_calc([zeros(length(v_trans_filtered),1),v_trans_filtered],0.15,0.10,60);
        else
            cls_onoff_det{c} = v_trans_filtered;
        end

        if contains(lock_where,'last')
            cls_onoff_det{c} = dur_block(ind_c) - (cls_onoff_det{c} - max(0,50-Md_temp));
        else
            cls_onoff_det{c} = cls_onoff_det{c} - max(0,50-Md_temp);
        end

        % if  contains(type_det , 'on')
        cls_onoff_det{c} = round(cls_onoff_det{c}) ; % diff correction
        % else
        %     cls_onoff_det{c} = round(cls_onoff_det{c}) ; % diff correction
        % end

        cls_onoff_det{c}(cls_onoff_det{c}<1) = 1;

    end


    temp = cell2mat(cls_onoff_det);
    temp(index_clustering(:,1)) = temp;
    cls_onoff_det_mat = [nan;temp;nan];

    this_fiducials = nan(1,length(ecg_rpeaks_index));
    for p = 2:length(ecg_rpeaks_index)-1
        this_index = bloks_index{p-1};
        this_fiducials(p) = this_index(min(length(this_index),cls_onoff_det_mat(p)));
        if p==2
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1)>0
                this_fiducials(1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1);
            end
        elseif p==length(ecg_rpeaks_index)-1
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1)<=length(data)
                this_fiducials(p+1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1);
            end
        end
    end

    positions.Poff = this_fiducials;


    lock_where = 'first';
    type_det = 'off';
    obs_seqment_beats = P_bloks_on;
    mean_init_1 = mean_init_Pon_1;
    mean_init_2 = mean_init_Pon_2;
    bloks_index = P_bloks_on_index;

    cls_onoff_det = cell(length(num_cls),1);
    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        cls_seqment_beats = obs_seqment_beats(:,ind_c); init_state_1 = mean_init_1(ind_c); init_state_2 = mean_init_2(ind_c);
        v_trans_det = lsim_fit_detect(cls_seqment_beats,type_det, init_state_1, init_state_2);

        if noisy_pwave(c)>0
            v_trans_det(:,1) = v_trans_det(:,2);
        end

        if contains(lock_where,'last')
            v_trans_det = dur_block(ind_c) - v_trans_det;
        end

        Md_temp = round(median(v_trans_det(:),'omitmissing'));
        v_trans_det = max(0,50-Md_temp) + v_trans_det;

        v_trans_m = movmedian(mean(v_trans_det,2),[100,100],'omitmissing');
        v_trans_det(:,3) = mean(v_trans_det,2);

        v_trans_filtered = v_trans_det(:,3);
        for p = 1:length(ind_c)
            if abs(v_trans_det(p,1)-v_trans_det(p,2))>sample_10ms
                [val_ms,ind_min] = min(abs(v_trans_det(p,:)-v_trans_m(p)));
                v_trans_filtered(p,1) = v_trans_det(p,ind_min);
            end
        end

        if flag_post_processing>0
            cls_onoff_det{c} = em_interval_calc([zeros(length(v_trans_filtered),1),v_trans_filtered],0.25,0.15,30);
        else
            cls_onoff_det{c} = v_trans_filtered;
        end

        if contains(lock_where,'last')
            cls_onoff_det{c} = dur_block(ind_c) - (cls_onoff_det{c} - max(0,50-Md_temp));
        else
            cls_onoff_det{c} = cls_onoff_det{c} - max(0,50-Md_temp);
        end

        % if  contains(type_det , 'on')
        cls_onoff_det{c} = round(cls_onoff_det{c}) ; % diff correction
        % else
        % cls_onoff_det{c} = round(cls_onoff_det{c}) ; % diff correction
        % end

        cls_onoff_det{c}(cls_onoff_det{c}<1) = 1;

    end


    temp = cell2mat(cls_onoff_det);
    temp(index_clustering(:,1)) = temp;
    cls_onoff_det_mat = [nan;temp;nan];

    this_fiducials = nan(1,length(ecg_rpeaks_index));
    for p = 2:length(ecg_rpeaks_index)-1
        this_index = bloks_index{p-1};
        this_fiducials(p) = this_index(min(length(this_index),cls_onoff_det_mat(p)));
        if p==2
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1)>0
                this_fiducials(1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(1);
            end
        elseif p==length(ecg_rpeaks_index)-1
            if this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1)<=length(data)
                this_fiducials(p+1) = this_fiducials(p)- ecg_rpeaks_index(p) + ecg_rpeaks_index(p+1);
            end
        end
    end

    positions.Pon = this_fiducials;
    positions.R = ecg_rpeaks_index(:)';

    % check the beat quality Pwave quality
    clear cls_mn_ecg
    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        cls_mn_ecg(c,:) = mean(ecg_blocks(ind_c,sample_70ms:end-sample_100ms),1);
    end

    beat_quality_score = (1/size(cls_mn_ecg,2)) *(ecg_blocks(:,sample_70ms:end-sample_100ms)-mean(ecg_blocks(:,sample_70ms:end-sample_100ms),2)) * (cls_mn_ecg-mean(cls_mn_ecg,2))' ./ (std(ecg_blocks(:,sample_70ms:end-sample_100ms),[],2)*std(cls_mn_ecg,[],2)');
    beat_quality_score = round(max(beat_quality_score,[],2),2);
    beat_quality_score = [0;beat_quality_score;0];

    QRSon_index = positions.QRSon;
    QRSon_index(isnan(QRSon_index)) = 1;

    QRS_dur = median(positions.QRSoff -  positions.QRSon,'omitmissing');
    R_amp = abs( median(data(positions.R)-data(QRSon_index)));
    P_thr = R_amp/20;

    P_indexes = [positions.Pon;positions.P; positions.Poff]';
    idx_notP = isnan(sum(P_indexes,2));
    P_indexes(idx_notP,:) = 1;

    P_amp = data(P_indexes(:,2))-(data(P_indexes(:,1))+data(P_indexes(:,3)))/2;
    sign_amp = sign(median(P_amp));

    if sign_amp>0
        P_amp = data(P_indexes(:,2))-max(data(P_indexes(:,1)),data(P_indexes(:,3)));
    else
        P_amp = data(P_indexes(:,2))-min(data(P_indexes(:,1)),data(P_indexes(:,3)));
    end

    P_amp_score = min(1,exp(3*sign_amp*(P_amp'-P_thr)/P_thr));
    P_dur_score = ( min(1,exp(-3*abs(P_indexes(:,3)-P_indexes(:,2)-QRS_dur/2)/max(sample_100ms/2,QRS_dur/2))) + min(1,exp(-3*abs(P_indexes(:,2)-P_indexes(:,1)-QRS_dur/2)/max(sample_100ms/2,QRS_dur/2))) )/2;
    PT_dur_score = (min(1,exp((min(sample_100ms, rr_intervals_ecg/6) -(QRSon_index'-P_indexes(:,1)))./min(sample_100ms, rr_intervals_ecg/6))) + min(1,exp(3*(QRSon_index'-P_indexes(:,3)-sample_10ms)./sample_10ms)) )/2;

    P_score = (P_amp_score + P_dur_score + PT_dur_score)/3;
    % P_score = movmean(P_score,[30,30]);
    P_score_fcm = P_score(:);
    P_score_fcm(isnan(P_score_fcm)) = 0.4;
    fcm_options = fcmOptions(NumClusters=2,MaxNumIteration=25, Exponent=1.1, Verbose=0);
    [fcm_centers, fcm_part_mat] = fcm(P_score_fcm(:),fcm_options);
    [~,cluster_fcm] = max(fcm_part_mat);
    [~,cls_max] = max(fcm_centers);
    thr_fcm = mean(P_score_fcm(cluster_fcm==cls_max)) - 3*std(P_score_fcm(cluster_fcm==cls_max));

    thr_fcm = max(thr_fcm,0.5);
    thr_fcm = min(thr_fcm,0.65);
    if flag_prune_P>0
        ind_dropP = find(P_score<thr_fcm); %find(movmin(P_score<0.5,[30,30]));
        positions.Pon(ind_dropP) = nan;
        positions.P(ind_dropP) = nan;
        positions.Poff(ind_dropP) = nan;
    end

    positions.P_score = P_score';
    positions.beat_quality_score = beat_quality_score';


    L = length(ecg_rpeaks_index_org);
    index_R = 1:L;
    index_R(index_remove) =[];
    positions_out.Pon = nan(1,L);    positions_out.Pon(index_R) = positions.Pon;
    positions_out.P = nan(1,L);      positions_out.P(index_R) = positions.P;
    positions_out.Poff = nan(1,L);   positions_out.Poff(index_R) = positions.Poff;

    positions_out.QRSon = nan(1,L);  positions_out.QRSon(index_R) = positions.QRSon;
    positions_out.R = nan(1,L);      positions_out.R = ecg_rpeaks_index_org(:)';
    positions_out.QRSoff = nan(1,L); positions_out.QRSoff(index_R) = positions.QRSoff;

    positions_out.Ton = nan(1,L);    positions_out.Ton(index_R) = positions.Ton;
    positions_out.T = nan(1,L);      positions_out.T(index_R) = positions.T;
    positions_out.Toff = nan(1,L);   positions_out.Toff(index_R) = positions.Toff;

    positions_out.P_score = zeros(1,L);  positions_out.P_score(index_R) = positions.P_score;
    positions_out.beat_quality_score = zeros(1,L); positions_out.beat_quality_score(index_R) = positions.beat_quality_score;

    positions = positions_out;

catch err

    warning('LSIM failed in computing all fiducial points')
    EXITFLAG.status = 'failed';
    EXITFLAG.message = err;
    L = length(ecg_rpeaks_index_org);
    index_R = 1:L;
    index_R(index_remove) =[];


    if ~isfield(positions,'Pon')
        positions_out.Pon = nan(1,L);
    else
        positions_out.Pon = nan(1,L);    positions_out.Pon(index_R) = positions.Pon;
    end
    if ~isfield(positions,'P')
        positions_out.P = nan(1,L);
    else
        positions_out.P = nan(1,L);      positions_out.P(index_R) = positions.P;
    end
    if ~isfield(positions,'Poff')
        positions_out.Poff = nan(1,L);
    else
        positions_out.Poff = nan(1,L);   positions_out.Poff(index_R) = positions.Poff;
    end

    if ~isfield(positions,'QRSon')
        positions_out.QRSon = nan(1,L);
    else
        positions_out.QRSon = nan(1,L); positions_out.QRSon(index_R) = positions.QRSon;
    end

    if ~isfield(positions,'R')
        positions_out.R = nan(1,L);
    else
        positions_out.R = nan(1,L);      positions_out.R = ecg_rpeaks_index_org(:)';
    end

    if ~isfield(positions,'QRSoff')
        positions_out.QRSoff = nan(1,L);
    else
        positions_out.QRSoff = nan(1,L); positions_out.QRSoff(index_R) = positions.QRSoff;
    end

    if ~isfield(positions,'Ton')
        positions_out.Ton = nan(1,L);
    else
        positions_out.Ton = nan(1,L);    positions_out.Ton(index_R) = positions.Ton;
    end

    if ~isfield(positions,'T')
        positions_out.T = nan(1,L);
    else
        positions_out.T = nan(1,L);      positions_out.T(index_R) = positions.T;
    end

    if ~isfield(positions,'Toff')
        positions_out.Toff = nan(1,L);
    else
        positions_out.Toff = nan(1,L);   positions_out.Toff(index_R) = positions.Toff;
    end


    if ~isfield(positions,'P_score')
        positions_out.P_score = zeros(1,L);
    else
        positions_out.P_score = zeros(1,L);  positions_out.P_score(index_R) = positions.P_score;
    end

    if ~isfield(positions,'beat_quality_score')
        positions_out.beat_quality_score = zeros(1,L);
    else
        positions_out.beat_quality_score = zeros(1,L); positions_out.beat_quality_score(index_R) = positions.beat_quality_score;
    end

    positions = positions_out;
end


end


function [beat_quality_score, ecg_rpeaks_index, P_rpeaks_index, max_min_r] = preprocess_rpeaks(data, ecg_rpeaks_index, fs, num_cluster)

fs_fd = fs;
sample_10ms = round(fs_fd*0.01);
sample_70ms = round(fs_fd*0.07);
sample_100ms = round(fs_fd*0.1);
sample_250ms = round(fs_fd*0.25);
sample_350ms = round(fs_fd*0.35);

rpeak_intervals = cell(length(ecg_rpeaks_index)-2,1);
for p = 2:length(ecg_rpeaks_index)-1
    rpeak_intervals{p-1} = (ecg_rpeaks_index(p)-sample_70ms:ecg_rpeaks_index(p)+sample_70ms)';
end
rpeak_intervals_mat = cell2mat( rpeak_intervals);
L = length(rpeak_intervals_mat);

for f = 1:8
    data_r = lp_filter_zero_phase(data, 35/fs);
    data_rh = data_r;
    data_r = data_r - lp_filter_zero_phase(data_r,f/fs);
    data_hat = [data_r(rpeak_intervals_mat)',ones(L,1)]*(pinv([data_r(rpeak_intervals_mat)',ones(L,1)])*data_rh(rpeak_intervals_mat)');

    temp_coef = corrcoef(data_hat,data_rh(rpeak_intervals_mat)');
    % Cxy(f) = 1-(norm(data_rh(rpeak_intervals_mat)' - data_hat)/  norm(data_rh(rpeak_intervals_mat)));
    Cxy(f) = temp_coef(1,2);
end
Cxy(f+1) =0;
f = find(Cxy<0.95,1,'first');
data_r = lp_filter_zero_phase(data, 35/fs);
data_r = data_r - lp_filter_zero_phase(data_r,f/fs);

P_rpeaks_index = ecg_rpeaks_index;
max_min_r = 0*ecg_rpeaks_index;

for p = 1:length(ecg_rpeaks_index)

    this_qrs_index = ecg_rpeaks_index(p)-sample_70ms:ecg_rpeaks_index(p)+sample_100ms;
    this_qrs_index(this_qrs_index<1) = [];
    this_qrs_index(this_qrs_index>length(data)) =[];

    [TF,P] = islocalmin(data_r(this_qrs_index),'MaxNumExtrema',1);
    TF_min = find(TF);
    if ~isempty(TF_min)
        P_min = P(TF_min);
    else
        P_min =0;
    end

    [TF,P] = islocalmax(data_r(this_qrs_index),'MaxNumExtrema',1);
    TF_max = find(TF);
    if ~isempty(TF_max)
        P_max = P(TF_max);
    else
        P_max =0;
    end


    if P_max>= P_min && P_max>0
        ecg_rpeaks_index(p) = this_qrs_index(TF_max);
        P_rpeaks_index(p) = P_max;
        max_min_r(p) = 1;
    elseif P_min>0
        ecg_rpeaks_index(p) = this_qrs_index(TF_min);
        P_rpeaks_index(p) = P_min;
        max_min_r(p) = -1;
    else

        [P_min,TF_min] = min(data_r(this_qrs_index));
        P_min = abs(P_min);

        [P_max,TF_max] = max(data_r(this_qrs_index));
        P_max = abs(P_max);



        if P_max>= P_min && P_max>0
            ecg_rpeaks_index(p) = this_qrs_index(TF_max);
            P_rpeaks_index(p) = P_max;
            max_min_r(p) = 1;
        elseif P_min>0
            ecg_rpeaks_index(p) = this_qrs_index(TF_min);
            P_rpeaks_index(p) = P_min;
            max_min_r(p) = -1;
        end


    end

end


avg_intervals_ecg = min(fs*1.5,1.25*median(diff(ecg_rpeaks_index),'omitmissing'));

pqrs_bloks = zeros(length(ecg_rpeaks_index)-2,sample_250ms+sample_70ms+1);
t_bloks = zeros(length(ecg_rpeaks_index)-2,round(sample_350ms*min(2,max(max(1,avg_intervals_ecg/(2*sample_350ms)))))-sample_70ms);
for p = 2:length(ecg_rpeaks_index)-1
    this_qrs_index = ecg_rpeaks_index(p)-min(sample_250ms,max(sample_70ms,floor(0.3*avg_intervals_ecg))):ecg_rpeaks_index(p)+sample_70ms;
    if any(this_qrs_index<1)||any(this_qrs_index>length(data))
        continue;
    end
    pqrs_bloks(p-1,end-length(this_qrs_index)+1:end) = data(this_qrs_index) ;


    this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:ecg_rpeaks_index(p)+min(sample_350ms,max(sample_250ms,floor(0.7*avg_intervals_ecg)));
    if any(this_t_index>length(data))
        continue;
    end

    t_bloks(p-1,1:length(this_t_index)) =  data(this_t_index);

end

ecg_blocks = [pqrs_bloks,t_bloks];

max_val = max(prctile(pqrs_bloks(:),99.9),prctile(t_bloks(:),99.9));
min_val = min(prctile(pqrs_bloks(:),0.1),prctile(t_bloks(:),0.1));

ecg_blocks(ecg_blocks>max_val) = max_val;
ecg_blocks(ecg_blocks<min_val) = min_val;
% ecg_blocks = (ecg_blocks - min_val) / (max_val-min_val);
for p = 1:size(ecg_blocks,1)
    ecg_blocks(p,:) = ecg_blocks(p,:) - linspace(mean(ecg_blocks(p,1:sample_10ms),2),mean(ecg_blocks(p,end-sample_10ms+1:end),2),size(ecg_blocks,2));
end

ecg_blocks(:,1:sample_10ms) = ecg_blocks(:,1:sample_10ms) .* linspace(0.01,1,sample_10ms);
ecg_blocks(:,end-sample_10ms+1:end) = ecg_blocks(:,end-sample_10ms+1:end) .* linspace(1,0.01,sample_10ms);

ecg_blocks_normalized = ecg_blocks./sqrt(sum(ecg_blocks.^2,2)); % normalization to unit power

fcm_features = ecg_blocks_normalized(:,sample_70ms:end);
fcm_features = fillmissing(fcm_features,"linear");
% fcm_features = zscore(fcm_features);

fcm_options = fcmOptions(NumClusters=num_cluster,MaxNumIteration=25, Exponent=1.1, Verbose=0);
[fcm_centers, fcm_part_mat] = fcm(fcm_features,fcm_options);
[~,cluster_fcm] = max(fcm_part_mat);

num_cls = size(fcm_centers,1);
index_clustering = [];
count_cls = zeros(num_cls,1);
for c = 1:num_cls
    index_temp = find(cluster_fcm==c);
    count_cls(c) = length(index_temp);
    index_clustering = [index_clustering;[index_temp(:),c*ones(length(index_temp),1)]];
end

if length(cluster_fcm)<30
    index_clustering(:,2)=1;
else
    index_clustering_org = index_clustering;
    for c = find(count_cls(:)'<20)
        [~,cls_max] = max(count_cls);
        index_clustering(index_clustering_org(:,2)==c,2) = cls_max;
    end
end

num_cls = unique(index_clustering(:,2));

% check the beat quality Pwave quality
clear cls_mn_ecg
for c = 1:length(num_cls)

    ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
    cls_mn_ecg(c,:) = mean(ecg_blocks(ind_c,sample_70ms:end-sample_100ms),1);
end

beat_quality_score = (1/size(cls_mn_ecg,2)) *(ecg_blocks(:,sample_70ms:end-sample_100ms)-mean(ecg_blocks(:,sample_70ms:end-sample_100ms),2)) * (cls_mn_ecg-mean(cls_mn_ecg,2))' ./ (std(ecg_blocks(:,sample_70ms:end-sample_100ms),[],2)*std(cls_mn_ecg,[],2)');
beat_quality_score = round(max(beat_quality_score,[],2),2);
beat_quality_score = [1;beat_quality_score;1];

end