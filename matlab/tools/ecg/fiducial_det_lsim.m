
function positions = fiducial_det_lsim(data, ecg_rpeaks_index, fs, varargin)

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
%        - win_T: Optional. It is window duration in second unit for extracting diff and std features from ecg signal for T-wave onset and offset (default is 0.02 or 20ms)
%        - win_P: Optional. It is window duration in second unit for extracting diff and std features from ecg signal for P-wave onset and offset (default is 0.02 or 20ms)
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
%   Reference:
%      ........
%
%   Sajjad Karimi, Reza Sameni  2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET
%
%


if nargin < 3
    error('The first 3 inputs are necessary for the fiducial detection')
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

fs_fd = 1000;
fs_fd = fs;

data = data(:)';
ecg_rpeaks_index = ecg_rpeaks_index(:);
ecg_rpeaks_index_org = ecg_rpeaks_index;
rr_intervals_ecg = em_interval_calc(ecg_rpeaks_index);
avg_intervals_ecg = movmean(rr_intervals_ecg,[60,60]);
index_remove = 1+find((diff(ecg_rpeaks_index(:))./avg_intervals_ecg-1)<-0.3 | diff(ecg_rpeaks_index(:))<0.25*fs);
ecg_rpeaks_index(index_remove) =[];

sample_10ms = round(fs_fd*0.01);
sample_70ms = round(fs_fd*0.07);
sample_100ms = round(fs_fd*0.1);
sample_250ms = round(fs_fd*0.25);
sample_350ms = round(fs_fd*0.35);

t_second = (0 : length(data)-1)/fs;

% Check optional input arguments
if nargin > 3 && ~isempty(varargin{1})
    flag_post_processing = varargin{1};
else
    flag_post_processing = 1;
end

if nargin > 4 && ~isempty(varargin{2})
    flag_prune_P = varargin{2};
else
    flag_prune_P = 1;
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

    fcm_features = [ecg_blocks_normalized(:,sample_70ms:end),repmat(avg_intervals_ecg(2:end-1),1,3*sample_10ms),repmat(linspace(-1,1,size(ecg_blocks,1))',1,3*sample_10ms)];
    fcm_features = fillmissing(fcm_features,"linear");
    fcm_features = zscore(fcm_features);

    fcm_options = fcmOptions(MaxNumIteration=25, Exponent=1.1, Verbose=0);
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


    ecg_denoised_std = movstd(ecg_denoised_n,[win_sample_qrs,win_sample_qrs]);
    ecg_denoised_ndiff = abs(ecg_denoised_ndiff);


    pqrs_bloks_on = cell(2,length(ecg_rpeaks_index)-2);
    pqrs_bloks_on_index = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_on_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_on_2 = cell(1,length(ecg_rpeaks_index)-2);

    for p = 2:length(ecg_rpeaks_index)-1

        max_sample_qrson = round((sample_100ms+2*sample_10ms)*max(0.9,rr_intervals_ecg(p-1)/(2*sample_350ms)));

        this_qrson_index = ecg_rpeaks_index(p)-max_sample_qrson:ecg_rpeaks_index(p);
        if any(this_qrson_index<1)
            pqrs_bloks_on{1,p-1} = ecg_denoised_ndiff(this_qrson_index);
            pqrs_bloks_on{1,p-1} = ecg_denoised_std(this_qrson_index) ;
            pqrs_bloks_on_index{p-1} = this_qrson_index;
            continue;
        end

        dur_this = round(length(this_qrson_index)/2);
        gain_sig = ones(1,length(this_qrson_index));
        ind_extrm =[];
        rpeak_val = abs(ecg_denoised_n(ecg_rpeaks_index(p)));
        if ecg_denoised_n(ecg_rpeaks_index(p))<0
            flag_type = 'max';
            [TF,P] = islocalmax(ecg_denoised_n(this_qrson_index(dur_this:end-sample_10ms)),'MaxNumExtrema',1,'MinProminence',0.1*rpeak_val);
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
            [TF,P] = islocalmin(ecg_denoised_n(this_qrson_index(dur_this:end-sample_10ms)),'MaxNumExtrema',1,'MinProminence',0.1*rpeak_val);
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
            temp = ecg_denoised_n(this_qrson_index(dur_this:end-sample_10ms));
            if contains(flag_type , 'max')
                sample_win_gain = TF - find(temp(1:TF)<temp(TF)-2*P/3,1,'last');
            else
                sample_win_gain = TF - find(temp(1:TF)>temp(TF)+2*P/3,1,'last');
            end
            temp = min(5,max(0,0.5*(rpeak_val/P)-1))*gausswin(2*sample_win_gain+1)';
            gain_sig(ind_extrm-sample_win_gain:ind_extrm) = gain_sig(ind_extrm-sample_win_gain:ind_extrm)+temp(sample_win_gain+1:end);
        end


        temp = ecg_denoised_std(this_qrson_index);
        thr_pq = mean(temp(end-3*sample_10ms:end))/5;
        [TF,P] = islocalmin(temp(1:end-2*sample_10ms),'MaxNumExtrema',1,'MinProminence',thr_pq/2);
        ind_drop_pq = find(TF);

        if isempty(ind_drop_pq)
            ind_drop_pq = 1;
        else
            ind_drop_pq = find(temp(1:ind_drop_pq)>thr_pq/2,1,'last');
        end

        this_qrson_index = this_qrson_index(ind_drop_pq:end);
        prior_win = gausswin(4*max(sample_70ms ,length(this_qrson_index))+1)';
        prior_win = prior_win(ceil(length(prior_win)/2)-2*length(this_qrson_index):ceil(length(prior_win)/2)+2*length(this_qrson_index));
        gain_sig = gain_sig(ind_drop_pq:end).*prior_win(length(this_qrson_index)+2:2*length(this_qrson_index)+1);

        pqrs_bloks_on{1,p-1} = ecg_denoised_ndiff(this_qrson_index).*gain_sig;
        pqrs_bloks_on{2,p-1} = ecg_denoised_std(this_qrson_index).*gain_sig ;

        pqrs_bloks_on_index{p-1} = this_qrson_index;
        dur_block(p-1,1) = length(this_qrson_index);

        ecg_denoised_ndiff(this_qrson_index) = ecg_denoised_ndiff(this_qrson_index).*gain_sig;
        ecg_denoised_std(this_qrson_index) = ecg_denoised_std(this_qrson_index).*gain_sig;

        dur_this = round(length(this_qrson_index)/2);
        mean_init_on_1{1,p-1} = [ecg_denoised_ndiff(this_qrson_index(1:dur_this));ecg_denoised_std(this_qrson_index(1:dur_this))];
        mean_init_on_2{1,p-1} = [ecg_denoised_ndiff(this_qrson_index(dur_this:end));ecg_denoised_std(this_qrson_index(dur_this:end))];

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
            cls_onoff_det{c} = round(cls_onoff_det{c})  - ceil(win_sample_qrs/2); % diff correction
        else
            cls_onoff_det{c} = round(cls_onoff_det{c}) + ceil(win_sample_qrs/2); % diff correction
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

    positions.QRSon = this_fiducials;


    ecg_qrson_index = positions.QRSon(:);
    pqrs_bloks_off = cell(2,length(ecg_rpeaks_index)-2);
    pqrs_bloks_off_index = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_off_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_off_2 = cell(1,length(ecg_rpeaks_index)-2);

    rcounter = 2;

    mn_qrson =  mean(ecg_rpeaks_index - ecg_qrson_index, 'omitmissing');
    sample_off_search = round(max(sample_100ms+2*sample_10ms, 2*sample_100ms - mn_qrson));

    for p = 2:length(ecg_rpeaks_index)-1

        this_qrsoff_index = ecg_rpeaks_index(p):ecg_rpeaks_index(p) + sample_off_search*min(2,max(1,avg_intervals_ecg(p)/(2*sample_350ms)));
        if any(this_qrsoff_index>length(data))
            continue;
        end

        dur_this = round(length(this_qrsoff_index)/2);
        gain_sig = ones(1,length(this_qrsoff_index));
        ind_extrm =[];
        rpeak_val = abs(ecg_denoised_n(ecg_rpeaks_index(p)));
        if ecg_denoised_n(ecg_rpeaks_index(p))<0
            flag_type = 'max';
            [TF,P] = islocalmax(ecg_denoised_n(ecg_qrson_index(p):ecg_rpeaks_index(p)),'MaxNumExtrema',1,'MinProminence',0.3*rpeak_val);
            TF = find(TF);
            if isempty(TF)
                rcounter = rcounter+1;
            end
            if isempty(TF) && rcounter/p>0.5

                [TF,P] = islocalmax(ecg_denoised_n(this_qrsoff_index(sample_10ms:dur_this)),'MaxNumExtrema',1,'MinProminence',0.1*rpeak_val);
                TF = find(TF);
                [p_sort,ind_p] = sort(P,"descend");
                if ~isempty(TF) || (p_sort(2)~=0 && p_sort(1)/p_sort(2)>10 && p_sort(1)>0.03*rpeak_val && ind_p(1)<ind_p(2)) || (p_sort(2)==0 && p_sort(1)>0.05*rpeak_val)
                    if isempty(TF)
                        [~,TF]= max(P);
                    end
                    ind_extrm = sample_10ms-1 +TF;
                    P = P(TF);
                end
            end
        else

            flag_type = 'min';
            [TF,P] = islocalmin(ecg_denoised_n(ecg_qrson_index(p):ecg_rpeaks_index(p)),'MaxNumExtrema',1,'MinProminence',0.3*rpeak_val);
            TF = find(TF);
            if isempty(TF)
                rcounter = rcounter+1;
            end
            if isempty(TF)&& rcounter/p>0.5

                [TF,P] = islocalmin(ecg_denoised_n(this_qrsoff_index(sample_10ms:dur_this)),'MaxNumExtrema',1,'MinProminence',0.1*rpeak_val);
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
        end

        if ~isempty(ind_extrm)
            temp = ecg_denoised_n(this_qrsoff_index(sample_10ms:dur_this));
            if contains(flag_type , 'max')
                sample_win_gain = find(temp(TF:end)<temp(TF)-2*P/3,1,'first');
            else
                sample_win_gain = find(temp(TF:end)>temp(TF)+2*P/3,1,'first');
            end
            temp = min(5,max(0,0.5*(rpeak_val/P)-1))*gausswin(2*sample_win_gain+1)';
            gain_sig(ind_extrm:ind_extrm+sample_win_gain) = gain_sig(ind_extrm:ind_extrm+sample_win_gain)+temp(1:sample_win_gain+1);
        end


        pqrs_bloks_off{1,p-1} = ecg_denoised_ndiff(this_qrsoff_index).*gain_sig;
        pqrs_bloks_off{2,p-1} = ecg_denoised_std(this_qrsoff_index).*gain_sig ;
        pqrs_bloks_off_index{p-1} = this_qrsoff_index;
        dur_block(p-1,1) = length(this_qrsoff_index);

        ecg_denoised_ndiff(this_qrsoff_index) = ecg_denoised_ndiff(this_qrsoff_index).*gain_sig;
        ecg_denoised_std(this_qrsoff_index) = ecg_denoised_std(this_qrsoff_index).*gain_sig;

        mean_init_off_1{1,p-1} = [ecg_denoised_ndiff(this_qrsoff_index(1:dur_this));ecg_denoised_std(this_qrsoff_index(1:dur_this))];
        mean_init_off_2{1,p-1} = [ecg_denoised_ndiff(this_qrsoff_index(dur_this:end));ecg_denoised_std(this_qrsoff_index(dur_this:end))];

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
            cls_onoff_det{c} = round(cls_onoff_det{c})  - ceil(win_sample_qrs/2); % diff correction
        else
            cls_onoff_det{c} = round(cls_onoff_det{c}) + ceil(win_sample_qrs/2); % diff correction
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

    % a = figure('Position', [130 130 1500 800]);
    % lg = {};
    % plot(t_second,ecg_denoised_n,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
    % hold on
    % plot(t_second,ecg_denoised_ndiff,LineWidth=1.5) ;lg = cat(2, lg, {'diff1'});
    % grid on
    % plot(t_second,ecg_denoised_std,LineWidth=1.5) ;lg = cat(2, lg, {'std'});

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
    ecg_denoised_nT = ecg_denoised_n;

    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_qrson_index(p)-min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1))) : ecg_qrsoff_index(p);
        ecg_denoised_nT(this_qrs_index) = linspace(ecg_denoised_nT(this_qrs_index(1)),ecg_denoised_nT(this_qrs_index(end)) , length(this_qrs_index) );

    end

    ecg_denoised_nT = lp_filter_zero_phase(ecg_denoised_nT, 10/fs);
    % ecg_denoised_nT = ecg_denoised_nT - lp_filter_zero_phase(ecg_denoised_nT, 0.5/fs);
    % ecg_denoised_nT = sjk_eeg_filter(ecg_denoised_nT, fs_fd,0.5,20);

    ecg_T_std = movstd(ecg_denoised_nT,[win_sample_T,win_sample_T]);
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
        temp(1)=ecg_denoised_nT(this_T_index(1));
        temp = temp  - linspace(temp(1),temp(end),length(temp));
        [TF,P] = islocalmax(temp,'MaxNumExtrema',2,'MinSeparation',sample_70ms);
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
        end

        TF_max = 5*sample_10ms+TF_max;


        [TF,P] = islocalmin(temp, 'MaxNumExtrema',2,'MinSeparation',sample_70ms);
        TF_min =find(TF);
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
            before_peak_cell{c,1} = em_interval_calc([zeros(length(ind_c),1), before_peak_cell{c,1}],0.15,0.10,60);
            before_peak_cell{c,1} = round( before_peak_cell{c,1})-max(0,50-Md_temp);

            after_peak_cell{c,1}  = after_peak(ind_c);
            Md_temp = round(median(after_peak_cell{c,1}(:),'omitmissing'));
            after_peak_cell{c,1} = max(0,50-Md_temp) + after_peak_cell{c,1};
            after_peak_cell{c,1} = em_interval_calc([zeros(length(ind_c),1),after_peak_cell{c,1}],0.15,0.10,60);
            after_peak_cell{c,1} = round(after_peak_cell{c,1})-max(0,50-Md_temp);

        end

        before_peak = cell2mat(before_peak_cell);
        before_peak(index_clustering(:,1)) = before_peak;
        before_peak = [nan;before_peak;nan];

        after_peak = cell2mat(after_peak_cell);
        after_peak(index_clustering(:,1)) = after_peak;
        after_peak = [nan;after_peak;nan];

    else
        before_peak = round([nan;before_peak;nan]);
        after_peak = round([nan;after_peak;nan]);
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


        % dur_this = ceil(length(this_bT_index)/3);

        % temp = ecg_T_std(this_bT_index);
        % if length(temp)>5*sample_10ms
        %     thr_tp = mean(temp(end-3*sample_10ms:end))/5;
        %     [TF,P] = islocalmin(temp(1:end-3*sample_10ms),'MaxNumExtrema',1,'MinProminence',thr_tp/2);
        %     ind_drop_pq = find(TF);
        %
        %     if isempty(ind_drop_pq)
        %         ind_drop_pq =  2* sample_10ms;
        %     else
        %         ind_drop_pq = ind_drop_pq - find(temp(1:ind_drop_pq)>thr_tp/2,1,'first');
        %     end
        %
        %     this_bT_index(1:ind_drop_pq)=[];
        % end
        % prior_win = gausswin(4*max(sample_70ms ,length(this_bT_index))+1)';
        % prior_win = prior_win(ceil(length(prior_win)/2)-2*length(this_bT_index):ceil(length(prior_win)/2)+2*length(this_bT_index));
        % gain_sig =prior_win(2*length(this_bT_index)+1:3*length(this_bT_index));

        T_bloks_on{1,p-1} = ecg_T_ndiff(this_bT_index);
        T_bloks_on{2,p-1} = ecg_T_std(this_bT_index);

        T_bloks_on_index{p-1} = this_bT_index;
        % ecg_T_ndiff(this_bT_index) = ecg_T_ndiff(this_bT_index).*gain_sig;
        % ecg_T_std(this_bT_index) = ecg_T_std(this_bT_index).*gain_sig;

        dur_this = ceil(length(this_bT_index)/3);

        mean_init_Ton_1{1,p-1} = [ecg_T_ndiff(this_bT_index(1:dur_this));ecg_T_std(this_bT_index(1:dur_this))];
        mean_init_Ton_2{1,p-1} = [ecg_T_ndiff(this_bT_index(dur_this:end));ecg_T_std(this_bT_index(dur_this:end))];

        this_aT_index = this_T_index(after_peak(p)+win_sample_T:end);

        if isempty(this_aT_index)
            this_aT_index = this_T_index(ceil(before_peak(p)/2):end);
        end


        temp = ecg_T_std(this_aT_index);
        if length(temp)>5*sample_10ms
            thr_tp = mean(temp(1:3*sample_10ms))/5;
            [TF,P] = islocalmin(temp(3*sample_10ms+1:end),'MaxNumExtrema',1,'MinProminence',thr_tp/2);
            ind_drop_pq = 3*sample_10ms+find(TF);

            if isempty(ind_drop_pq)
                ind_drop_pq = length(temp);
            else
                ind_drop_pq = ind_drop_pq + find(temp(ind_drop_pq+1:end)>thr_tp/2,1,'first');
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
            cls_onoff_det{c} = round(cls_onoff_det{c}) - ceil(win_sample_T/2); % diff correction
        else
            cls_onoff_det{c} = round(cls_onoff_det{c}); % diff correction
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
            cls_onoff_det{c} = round(cls_onoff_det{c}) - ceil(win_sample_T/2); % diff correction
        else
            cls_onoff_det{c} = round(cls_onoff_det{c}); % diff correction
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


    % a = figure('Position', [130 130 1500 800]);
    % lg = {};
    % plot(t_second,ecg_denoised_n,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
    % hold on
    % plot(t_second,ecg_denoised_nT,LineWidth=1.5) ;lg = cat(2, lg, {'ECG-T'});
    % plot(t_second,ecg_T_ndiff,LineWidth=1.5) ;lg = cat(2, lg, {'diff1'});
    % % plot(t_second,ecg_T_ndiff2,LineWidth=1.5) ;lg = cat(2, lg, {'diff2'});
    % grid on
    % plot(t_second,ecg_T_std,LineWidth=1.5) ;lg = cat(2, lg, {'std'});

    %% P-wave Detection

    P_bloks_off = cell(2,length(ecg_rpeaks_index)-2);
    P_bloks_off_index = cell(1,length(ecg_rpeaks_index)-2);

    P_bloks_on = cell(2,length(ecg_rpeaks_index)-2);
    P_bloks_on_index = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_Poff_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_Poff_2 = cell(1,length(ecg_rpeaks_index)-2);

    mean_init_Pon_1 = cell(1,length(ecg_rpeaks_index)-2);
    mean_init_Pon_2 = cell(1,length(ecg_rpeaks_index)-2);

    ecg_Toff_index = positions.Toff;
    ecg_denoised_nP = 0*data;

    max_indexes = nan(length(ecg_rpeaks_index)-2,1);
    max_P_vals = nan(length(ecg_rpeaks_index)-2,1);
    min_indexes = nan(length(ecg_rpeaks_index)-2,1);
    min_P_vals = nan(length(ecg_rpeaks_index)-2,1);
    P_wave_blocks = zeros(length(ecg_rpeaks_index)-2, 2*sample_250ms);
    for p = 2:length(ecg_rpeaks_index)-1

        temp_back = min( 2*sample_250ms, max(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
        this_P_index = ecg_qrson_index(p)-temp_back: ecg_qrson_index(p)-sample_10ms;
        if any(this_P_index<1)||any(this_P_index>length(data))
            continue;
        end

        ecg_denoised_nP(this_P_index) = data(this_P_index) - linspace(data(this_P_index(1)),data(this_P_index(end)) , length(this_P_index) );
        P_wave_blocks(p-1,end-length(this_P_index)+1:end) = ecg_denoised_nP(this_P_index);
        ecg_denoised_nP(this_P_index) = ecg_denoised_nP(this_P_index)/norm(ecg_denoised_nP(this_P_index));

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

    if abs(mean(max_P_vals,'omitmissing')) > abs(mean(min_P_vals,'omitmissing'))
        std_index = max_indexes;
    else
        std_index = min_indexes;
    end

    noisy_pwave = 0;
    clear mn_std_index

    for c = 1:length(num_cls)

        ind_c = index_clustering(index_clustering(:,2)==num_cls(c),1);
        mn_std_index(c) = mean(movstd(std_index(ind_c),[60,60]),'omitmissing');

        if mean(mn_std_index(c),'omitmissing')>sample_10ms
            noisy_pwave(c) = 1;
            P_wave_blocks(ind_c,:) = movmean(P_wave_blocks(ind_c,:),[30,30]);
        else
            noisy_pwave(c) = 0;
        end

    end

    for p = 2:length(ecg_rpeaks_index)-1

        temp_back = min( 2*sample_250ms, max(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
        this_P_index = ecg_qrson_index(p)-temp_back: ecg_qrson_index(p)-sample_10ms;
        if any(this_P_index<1)||any(this_P_index>length(data))
            continue;
        end

        ecg_denoised_nP(this_P_index) = P_wave_blocks(p-1,end-length(this_P_index)+1:end) ;
        ecg_denoised_nP(this_P_index) = ecg_denoised_nP(this_P_index)/norm(ecg_denoised_nP(this_P_index));
    end

    ecg_P_std = movstd(ecg_denoised_nP,[ceil(win_sample_P/2),ceil(win_sample_P/2)]);
    ecg_P_ndiff = [zeros(1,win_sample_P-ceil(win_sample_P/2)), ecg_denoised_nP(win_sample_P+1:end) - ecg_denoised_nP(1:end-win_sample_P),zeros(1,win_sample_P-floor(win_sample_P/2))];

    ecg_P_ndiff = abs(ecg_P_ndiff);

    % a = figure('Position', [130 130 1500 800]);
    % lg = {};
    % plot(t_second,ecg_denoised_n,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
    % hold on
    % plot(t_second,ecg_denoised_nP,LineWidth=1.5) ;lg = cat(2, lg, {'ECG-T'});
    % plot(t_second,ecg_P_ndiff,LineWidth=1.5) ;lg = cat(2, lg, {'diff1'});
    % grid on
    % plot(t_second,ecg_P_std,LineWidth=1.5) ;lg = cat(2, lg, {'std'});


    P_peaks = nan(length(ecg_rpeaks_index),1);
    Pbefore_peak = nan(length(ecg_rpeaks_index),1);
    Pafter_peak = nan(length(ecg_rpeaks_index),1);
    Pbiphasic_shape = nan(length(ecg_rpeaks_index),1);
    P_type = cell(length(ecg_rpeaks_index),1);
    P_type{1} = 'x';
    P_type{end} = 'x';

    for p = 2:length(ecg_rpeaks_index)-1

        % this_P_index = ecg_qrsoff_index(p)-round(min(0.66*(ecg_Toff_index(p-1)-ecg_qrsoff_index(p))),  (sample_100ms+sample_70ms)*max(1,avg_intervals_ecg(p)/(2*sample_350ms))):ecg_qrsoff_index(p);

        temp_back = min( 2*sample_250ms, max(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
        this_P_index = ecg_qrson_index(p)-temp_back: ecg_qrson_index(p);

        if any(this_P_index>length(data))
            continue;
        end

        if length(this_P_index)<2*sample_10ms
            this_P_index = ecg_qrson_index(p)-sample_100ms-sample_70ms: ecg_qrson_index(p);
        end

        this_P_index = flip(this_P_index);
        [TF,P] = islocalmax(ecg_denoised_nP(this_P_index),'MaxNumExtrema',2,'MinSeparation',sample_10ms);
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
        end

        [TF,P] = islocalmin(ecg_denoised_nP(this_P_index),'MaxNumExtrema',2,'MinSeparation',sample_10ms);
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

    Md_temp = round(median(Pbefore_peak(:),'omitmissing'));
    Pbefore_peak = max(0,50-Md_temp) + Pbefore_peak;
    Pbefore_peak = em_interval_calc([zeros(length(ecg_rpeaks_index),1),Pbefore_peak],0.15,0.10,60);
    Pbefore_peak = round(Pbefore_peak)-max(0,50-Md_temp);

    Md_temp = round(median(Pafter_peak(:),'omitmissing'));
    Pafter_peak = max(0,50-Md_temp) + Pafter_peak;
    Pafter_peak = em_interval_calc([zeros(length(ecg_rpeaks_index),1),Pafter_peak],0.15,0.10,60);
    Pafter_peak = round(Pafter_peak)-max(0,50-Md_temp);

    for p = 2:length(ecg_rpeaks_index)-1

        if isnan(Pafter_peak(p))
            Pafter_peak(p) = 1;
        end

        temp_back = min( 2*sample_250ms, max(ecg_qrson_index(p)-ecg_Toff_index(p-1)-sample_10ms ,min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1)))));
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
        elseif contains(P_type{p},'bi-phasic')
            P_peaks(p) = this_P_index(Pbefore_peak(p));
        elseif contains(P_type{p},'Umax') ||  contains(P_type{p},'Umin')
            P_peaks(p) = this_P_index(Pbefore_peak(p));
            this_P_index = this_P_index(1: min(T_dur,max(Pbefore_peak(p)+3*sample_10ms,Pbiphasic_shape(p))));
        elseif contains(P_type{p},'Nmax') ||  contains(P_type{p},'Nmin')
            P_peaks(p) = this_P_index(Pbefore_peak(p));
            this_P_index = this_P_index(1: min(T_dur,Pbiphasic_shape(p)));
        end

        try
            this_bP_index = this_P_index(1:Pbefore_peak(p));
        catch
            this_bP_index = this_P_index(1:round(Pbefore_peak(p)/2));
        end

        % dur_this = ceil(length(this_bP_index)/3);
        % thr_PQ = max(ecg_P_std(this_bP_index(2*dur_this:end)))/2;
        % ind_drop_st = min(dur_this,find(ecg_P_std(this_bP_index)<thr_PQ,1,"first"));
        % if isempty(ind_drop_st)
        %     ind_drop_st = sample_10ms;
        % end
        % ratio_std = ecg_P_std(this_bP_index(ind_drop_st))/thr_PQ;
        % this_bP_index = this_bP_index(ind_drop_st:end);
        % if ratio_std<0.5
        %     P_bloks_off{1,p-1} = ecg_P_ndiff(this_bP_index);
        %     P_bloks_off{2,p-1} = ecg_P_std(this_bP_index);
        % else
        %     P_bloks_off{1,p-1} = ecg_P_ndiff(this_bP_index); P_bloks_off{1,p-1}(:,1:sample_10ms) = P_bloks_off{1,p-1}(:,1:sample_10ms).*rand(1,sample_10ms);
        %     P_bloks_off{2,p-1} = ecg_P_std(this_bP_index);  P_bloks_off{2,p-1}(:,1:sample_10ms) = P_bloks_off{2,p-1}(:,1:sample_10ms).*rand(1,sample_10ms);
        % end

        P_bloks_off{1,p-1} = ecg_P_ndiff(this_bP_index); %P_bloks_on{1,p-1}(:,1:2*sample_10ms) = P_bloks_on{1,p-1}(:,1:2*sample_10ms).*linspace(0,1,2*sample_10ms);
        P_bloks_off{2,p-1} = ecg_P_std(this_bP_index);  %P_bloks_on{2,p-1}(:,1:2*sample_10ms) = P_bloks_on{2,p-1}(:,1:2*sample_10ms).*linspace(0,1,2*sample_10ms);
        P_bloks_off_index{p-1} = this_bP_index;

        dur_this = ceil(length(this_bP_index)/3);
        mean_init_Poff_1{1,p-1} = [ecg_P_ndiff(this_bP_index(1:dur_this));ecg_P_std(this_bP_index(1:dur_this))];
        mean_init_Poff_2{1,p-1} = [ecg_P_ndiff(this_bP_index(dur_this:end));ecg_P_std(this_bP_index(dur_this:end))];

        this_aP_index = this_P_index(Pafter_peak(p):end);

        if isempty(this_aP_index)
            this_aP_index = this_P_index(ceil(Pbefore_peak(p)/2):end);
        end
        P_bloks_on{1,p-1} = ecg_P_ndiff(this_aP_index);
        P_bloks_on{2,p-1} = ecg_P_std(this_aP_index);
        P_bloks_on_index{p-1} = this_aP_index;

        dur_this = ceil(length(this_aP_index)/2);
        mean_init_Pon_1{1,p-1} = [ecg_P_ndiff(this_aP_index(1:dur_this));ecg_P_std(this_aP_index(1:dur_this))];
        mean_init_Pon_2{1,p-1} = [ecg_P_ndiff(this_aP_index(dur_this:end));ecg_P_std(this_aP_index(dur_this:end))];


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
    P_score = movmean(P_score,[30,30]);

    if flag_prune_P>0
        ind_dropP = find(movmin(P_score<0.5,[30,30]));
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

catch


    L = length(ecg_rpeaks_index_org);

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

end


end

