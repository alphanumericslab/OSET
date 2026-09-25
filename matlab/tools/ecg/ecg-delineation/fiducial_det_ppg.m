
function [positions, EXITFLAG] = fiducial_det_ppg(data, ecg_rpeaks_index, fs, varargin)

% ECG fiducial points detector based on LSIM approach
%   [positions, EXITFLAG] = fiducial_det_lsim(data, ecg_rpeaks_index, fs, varargin)
%
%   Inputs:
%       data: Vector of input ECG data
%       ecg_rpeaks_index: R-peak indexes in samples
%       fs: Sampling rate in Hz
%  (optional input):
%        - flag_post_processing: Optional. It is a flag that can be 0 or 1, if it is set to 1 then prioir
%                   information from RR intervals is involved to modifying initail estimated fiducial points by LSIM (default is 0)
%        - flag_prune_P: Optional. It is a flag that can be 0 or 1, if it is set
%                   to 1 then based on a score for p-wave quality assesment the detected p-waves are pruned  (default is 0)
%        - win_qrs: Optional. It is window duration in second unit for extracting diff and std features from ecg signal for QRS onset and offset (default is 0.01 or 10ms)
%        - win_T: Optional. It is window duration in second unit for extracting diff and std features from ecg signal for T-wave onset and offset (default is 0.02 or 20ms)
%        - win_P: Optional. It is window duration in second unit for extracting diff and std features from ecg signal for P-wave onset and offset (default is 0.02 or 20ms)
%        - num_cluster: Number of clusters for FCM (default is 3)
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


fs_fd = fs;

sample_10ms = round(fs_fd*0.01);
sample_70ms = round(fs_fd*0.07);
sample_350ms = round(fs_fd*0.35);
sample_550ms = round(fs_fd*0.55);

% Check optional input arguments
if nargin > 3 && ~isempty(varargin{1})
    flag_post_processing = varargin{1};
else
    flag_post_processing = 0;
end


if length(ecg_rpeaks_index)<3

    num_beats = length(ecg_rpeaks_index);
    positions = struct('R',nan(1,num_beats),'onset',nan(1,num_beats),'sys_peak',nan(1,num_beats),...
        'd_notch',nan(1,num_beats),'dias_peak',nan(1,num_beats),'u_peak',nan(1,num_beats),...
        'w_peak',nan(1,num_beats));

    return
end


data = data(:)';

for f_hp= 0.5:0.1:1.5

    ppg_raw = data(:)';
    ppg_raw = ppg_raw - (movmean(movmedian(ppg_raw(:),[floor(2*fs),floor(2*fs)]),[floor(1*fs),floor(1*fs)]))';
    ppg_denoised = lp_filter_zero_phase(ppg_raw, 2*8/fs);
    ppg_denoised = ppg_denoised - lp_filter_zero_phase(ppg_denoised, 2*f_hp/fs);
    ppg_denoised = ppg_denoised./(movstd(ppg_denoised(:),[floor(5*fs),floor(5*fs)]))';

    [a,b]=pwelch(ppg_denoised(:),round(4*fs),round(2*fs),round(4*fs),fs);
    ind_15 = find(b<=2,1,'last');
    ind_10 = find(b<=5,1,'last');
    if sum(a(1:ind_15(1)))/sum(a(1:ind_10(1)))<0.25
        % disp(sum(a(1:ind_15))/sum(a(1:ind_10)))
        % disp(['fhp: ',num2str(f_hp)])
        break;
    end

end

% f_hp = 2;
positions_out.f_hp = f_hp;

ppg_raw = data(:)';
ppg_raw = ppg_raw - (movmean(movmedian(ppg_raw(:),[floor(2*fs),floor(2*fs)]),[floor(1*fs),floor(1*fs)]))';
ppg_denoised = lp_filter_zero_phase(ppg_raw, 2*8/fs);
ppg_denoised = ppg_denoised - lp_filter_zero_phase(ppg_denoised, 2*f_hp/fs);
data = ppg_denoised./(movstd(ppg_denoised(:),[floor(5*fs),floor(5*fs)]))';
dppg = [0, (diff(data(:)))'];
dppg = 0.5*dppg./(movstd(dppg(:),[floor(5*fs),floor(5*fs)]))';


data_cwd = data - lp_filter_zero_phase(data, 2*(f_hp+0.25)/fs);
dppg_cwd = [0, (diff(data_cwd(:)))'];
dppg_cwd = 0.5*dppg_cwd./(movstd(dppg_cwd(:),[floor(5*fs),floor(5*fs)]))';

ecg_rpeaks_index = ecg_rpeaks_index(:);
ecg_rpeaks_index_org = ecg_rpeaks_index;
ind_org = 1:length(ecg_rpeaks_index_org);


rr_intervals_ecg = em_interval_calc(ecg_rpeaks_index);
N = max(5,min(60,ceil(length(ecg_rpeaks_index)/20)));
rr_intervals_ecg(rr_intervals_ecg>1.75*fs) = 1.75*fs;
avg_intervals_ecg = movmean(rr_intervals_ecg,[N,N]);
index_remove = 1+find(((diff(ecg_rpeaks_index(:))./avg_intervals_ecg-1)<-0.4 & ~(diff(ecg_rpeaks_index(:))>0.8*fs)) | diff(ecg_rpeaks_index(:))<0.25*fs);
ecg_rpeaks_index(index_remove) =[];

index_remove2 = ind_org(index_remove);
index_remove = unique([index_remove2(:)]);

num_beats = length(ecg_rpeaks_index);
positions = struct('R',nan(1,num_beats),'onset',nan(1,num_beats),'sys_peak',nan(1,num_beats),...
    'd_notch',nan(1,num_beats),'dias_peak',nan(1,num_beats),'u_peak',nan(1,num_beats),...
    'w_peak',nan(1,num_beats));


% G = gcd( fs_fd,fs );
% data = resample(data,fs_fd/G,fs/G);

% try

% rr_intervals_ecg = em_interval_calc(ecg_rpeaks_index);
% N = min(60,ceil(num_beats/20));
% rr_intervals_ecg(rr_intervals_ecg>1.5*fs) = 1.5*fs;
% avg_intervals_ecg = movmean(rr_intervals_ecg,[N,N]);
% avg_intervals_ecg = [avg_intervals_ecg;avg_intervals_ecg(end)];
% rr_intervals_ecg = [rr_intervals_ecg;rr_intervals_ecg(end)];

ecg_rpeaks_index = ecg_rpeaks_index(:)';
sys_peak_index = nan*ecg_rpeaks_index;
sys_peak_index_max = nan*ecg_rpeaks_index;
du_peak_index = nan*ecg_rpeaks_index;
dw_peak_index = nan*ecg_rpeaks_index;
onset_index = nan*ecg_rpeaks_index;

dnotch_peak_index = nan*ecg_rpeaks_index;
dias_peak_index = nan*ecg_rpeaks_index;

thr_5perc = prctile(data,5);
for p = 1:num_beats-1

    % finding Systolic peak
    temp_index = ecg_rpeaks_index(p)+sample_70ms:min(ecg_rpeaks_index(p+1)-sample_70ms,ecg_rpeaks_index(p)+5*sample_350ms);

    [TF,P_sys] = islocalmax(data(temp_index),'MaxNumExtrema',1);
    TF = find(TF);
    if isempty(TF)
        continue;
    end
    P_sys = P_sys(TF(1));
    sys_peak_index(p) = temp_index(TF(1));
    sys_peak_index_max(p) = temp_index(TF(1));

    % finding U peak in first order diff
    temp_index = ecg_rpeaks_index(p)+sample_10ms:sys_peak_index_max(p);
    [TF,~] = islocalmax(dppg(temp_index),'MaxNumExtrema',1);
    TF = find(TF);
    if isempty(TF)
        continue;
    end
    du_peak_index(p) = temp_index(TF(1));

    % finding onset of PPG
    temp_index = ecg_rpeaks_index(p)+sample_10ms: du_peak_index(p);
    TF = find(dppg(temp_index).*linspace(0.05,1,length(temp_index))<0.1*dppg(du_peak_index(p)),1,'last');
    if isempty(TF)
        continue;
    end
    onset_index(p) = temp_index(TF(1));

    % Correcting Systolic peak as first local peak
    % temp_index = sys_peak_index_max(p)+sample_10ms:sys_peak_index_max(p)+2*sample_70ms;
    % [TF,P_sys2] = islocalmax(data(temp_index),'MaxNumExtrema',1);
    % TF = find(TF);

    temp_index = du_peak_index(p):sys_peak_index_max(p)-2*sample_10ms;
    TF = find(dppg(temp_index)<0.01*dppg(du_peak_index(p)),1,'first');
    % TF2 = find(dppg(temp_index)<0.05*dppg(du_peak_index(p)),1,'first');
    if ~isempty(TF) && (sys_peak_index_max(p) -onset_index(p) > 2*(sample_70ms+sample_10ms)) % && (data(temp_index(TF2(1)))>0.9*data(sys_peak_index_max(p))||data(temp_index(TF(1)))>0.8*data(sys_peak_index_max(p)))
        sys_peak_index(p) = temp_index(TF(1));
    else
        temp_index = du_peak_index(p):sys_peak_index_max(p)+sample_10ms;
        [TF,P_n] = islocalmin(dppg(temp_index),'MaxNumExtrema',1);
        TF = find(TF);
        if ~isempty(TF) && P_n(TF)>0.05*(dppg(du_peak_index(p)) - dppg(temp_index(TF(1))))
            sys_peak_index(p) = temp_index(TF(1));
        end
        % sys_peak_index_max(p) = temp_index(TF(1));
    end

    % finding Dicrotic Notch
    temp_index = max(onset_index(p)+3*sample_70ms, sys_peak_index_max(p)+sample_70ms):min( onset_index(p)+sample_550ms , ecg_rpeaks_index(p+1)-3*sample_10ms);
    [TFdn,Pdn] = islocalmin(data(temp_index),'MaxNumExtrema',1,'MinProminence',P_sys/50);
    TFdn = find(TFdn);

    [TFda,Pda] = islocalmax(data(temp_index),'MaxNumExtrema',1,'MinProminence',P_sys/50);
    TFda = find(TFda);


    if ~isempty(TFdn) && ~isempty(TFda) && TFda< TFdn && Pdn(TFdn)<0.75*Pda(TFda)
        temp_index = max(onset_index(p)+2*(sample_70ms+sample_10ms), sys_peak_index_max(p)+sample_70ms-2*sample_10ms):min( onset_index(p)+sample_550ms , ecg_rpeaks_index(p+1)-3*sample_10ms);
        [TFdn,~] = islocalmin(data(temp_index),'MaxNumExtrema',1,'MinProminence',P_sys/50);
        TFdn = find(TFdn);
    end

    if isempty(TFdn) || data(temp_index(TFdn))<thr_5perc
        [TFdw,Pdw] = islocalmax(dppg_cwd(temp_index),'MaxNumExtrema',2);
        TFdw = find(TFdw);
        temp = dppg_cwd(temp_index);
        if ~isempty(TFdw)
            if ~isempty(TFdw)
                P = Pdw;
                if isscalar(TFdw)
                    % P_min= P(TFdw) + temp(TFdw);
                    P(TFdw) = 0;
                    [~, TF_minU] = max(P);
                    TFdw(2) = TF_minU;
                else
                    P_min=P(TFdw(1)) + temp(TFdw(1));
                    P_minU=P(TFdw(2))+ temp(TFdw(2));
                    if P_min < P_minU && P(TFdw(1))/P(TFdw(2))<2
                        TFdw = flip(TFdw);
                        % P_min=P(TFdw(1))+ abs(temp(TFdw(1)));
                        % P_minU=P(TFdw(2))+ abs(temp(TFdw(2)));
                    end
                end
            end
            dw_peak_index(p) = temp_index(TFdw(1));
            temp_index = max(sys_peak_index_max(p),dw_peak_index(p)-2*sample_70ms):min(ecg_rpeaks_index(p+1)-sample_10ms,dw_peak_index(p)+2*sample_70ms);
            % data_detrended =  data_cwd(temp_index);
            data_detrended =  data_cwd(temp_index) - linspace( data_cwd(temp_index(1)),data_cwd(temp_index(end)), length(temp_index));
            [TFdn,~] = islocalmin(data_detrended,'MaxNumExtrema',1);
            TFdn = find(TFdn);

            [TFdias,~] = islocalmax(data_detrended,'MaxNumExtrema',1);
            TFdias = find(TFdias);

            if ~isempty(TFdn) &&  ~isempty(TFdias) && TFdias > TFdn
                dnotch_peak_index(p) = temp_index(TFdn(1));
                dias_peak_index(p) = temp_index(TFdias);
            end

        end

    else
        dnotch_peak_index(p) = temp_index(TFdn(1));
        temp_index = temp_index(max(1,TFdn(1)-2*sample_10ms):min(length(temp_index),max(TFdn(1)+2*sample_70ms+3*sample_10ms)));
        [TFdw,~] = islocalmax(dppg(temp_index),'MaxNumExtrema',1);
        TFdw = find(TFdw);
        if ~isempty(TFdw)
            dw_peak_index(p) = temp_index(TFdw(1));
            temp_index = dw_peak_index(p): min(ecg_rpeaks_index(p+1)-sample_10ms, dw_peak_index(p)+ 3*(dw_peak_index(p)-dnotch_peak_index(p)));
            [TFdias,~] = islocalmax(data_cwd(temp_index),'MaxNumExtrema',1);
            TFdias = find(TFdias);
            if ~isempty(TFdias)
                dias_peak_index(p) = temp_index(TFdias);
            end
        else
            dnotch_peak_index(p) = nan;
        end
    end


end



positions.onset = onset_index;
positions.sys_peak = sys_peak_index;
positions.d_notch = dnotch_peak_index;
positions.dias_peak = dias_peak_index;
positions.u_peak = du_peak_index;
positions.w_peak = dw_peak_index;

if flag_post_processing>0

    % Refine onset position relative to systolic peak
    this_fiducials = positions.onset;
    dur_qrs = sys_peak_index(:) - this_fiducials(:);
    md_temp = median(dur_qrs);
    dur_qrs = max(0,0.1*fs-md_temp) + dur_qrs;

    dur_qrs = em_interval_calc([zeros(length(dur_qrs),1),dur_qrs],0.20,0.15,180);
    ind_median = abs(dur_qrs - movmedian(dur_qrs,[120,120]))/median(dur_qrs)>0.1;
    dur_qrs(ind_median) = median(dur_qrs);

    dur_qrs = round(dur_qrs-max(0,0.1*fs-md_temp));
    positions.onset = sys_peak_index(:)' - dur_qrs(:)'; % Onset position before systolic peak

    % Refine u_peak position relative to systolic peak
    this_fiducials = positions.u_peak;
    dur_qrs = sys_peak_index(:) - this_fiducials(:);
    md_temp = median(dur_qrs);
    dur_qrs = max(0,0.1*fs-md_temp) + dur_qrs;

    dur_qrs = em_interval_calc([zeros(length(dur_qrs),1),dur_qrs],0.20,0.15,180);
    ind_median = abs(dur_qrs - movmedian(dur_qrs,[120,120]))/median(dur_qrs)>0.1;
    dur_qrs(ind_median) = median(dur_qrs);

    dur_qrs = round(dur_qrs-max(0,0.1*fs-md_temp));
    positions.u_peak = sys_peak_index(:)' - dur_qrs(:)'; % U wave peak position before systolic peak

    % Refine w_peak position relative to systolic peak
    this_fiducials = positions.w_peak;
    dur_qrs = this_fiducials(:) - sys_peak_index(:);
    md_temp = median(dur_qrs);
    dur_qrs = max(0,0.1*fs-md_temp) + dur_qrs;

    dur_qrs = em_interval_calc([zeros(length(dur_qrs),1),dur_qrs],0.20,0.15,180);
    ind_median = abs(dur_qrs - movmedian(dur_qrs,[120,120]))/median(dur_qrs)>0.1;
    dur_qrs(ind_median) = median(dur_qrs);

    dur_qrs = round(dur_qrs-max(0,0.1*fs-md_temp));
    dur_qrs = round(dur_qrs-max(0,0.1*fs-md_temp));
    positions.w_peak  = sys_peak_index(:)' + dur_qrs(:)'; % W wave peak position after systolic peak

    % Refine dicrotic notch position relative to systolic peak
    this_fiducials = positions.d_notch;
    dur_qrs = this_fiducials(:) - sys_peak_index(:);
    md_temp = median(dur_qrs);
    dur_qrs = max(0,0.1*fs-md_temp) + dur_qrs;

    dur_qrs = em_interval_calc([zeros(length(dur_qrs),1),dur_qrs],0.20,0.15,180);
    ind_median = abs(dur_qrs - movmedian(dur_qrs,[120,120]))/median(dur_qrs)>0.1;
    dur_qrs(ind_median) = median(dur_qrs);

    dur_qrs = round(dur_qrs-max(0,0.1*fs-md_temp));
    dur_qrs = round(dur_qrs-max(0,0.1*fs-md_temp));
    positions.d_notch  = sys_peak_index(:)' + dur_qrs(:)'; % Dicrotic notch position after systolic peak

    % Refine diastolic peak position relative to systolic peak
    this_fiducials = positions.dias_peak;
    dur_qrs = this_fiducials(:) - sys_peak_index(:);
    md_temp = median(dur_qrs);
    dur_qrs = max(0,0.1*fs-md_temp) + dur_qrs;

    dur_qrs = em_interval_calc([zeros(length(dur_qrs),1),dur_qrs],0.20,0.15,180);
    ind_median = abs(dur_qrs - movmedian(dur_qrs,[120,120]))/median(dur_qrs)>0.1;
    dur_qrs(ind_median) = median(dur_qrs);

    dur_qrs = round(dur_qrs-max(0,0.1*fs-md_temp));
    dur_qrs = round(dur_qrs-max(0,0.1*fs-md_temp));
    positions.dias_peak  = sys_peak_index(:)' + dur_qrs(:)'; % Diastolic peak position after systolic peak


end

positions.R = ecg_rpeaks_index;

L = length(ecg_rpeaks_index_org);
index_R = 1:L;
index_R(index_remove) =[];

positions_out.R = nan(1,L);      positions_out.R = ecg_rpeaks_index_org(:)';
positions_out.onset = nan(1,L);    positions_out.onset(index_R) = positions.onset;
positions_out.u_peak = nan(1,L);   positions_out.u_peak(index_R) = positions.u_peak;
positions_out.sys_peak = nan(1,L);      positions_out.sys_peak(index_R) = positions.sys_peak;

positions_out.d_notch = nan(1,L);      positions_out.d_notch(index_R) = positions.d_notch;
positions_out.w_peak = nan(1,L);  positions_out.w_peak(index_R) = positions.w_peak;
positions_out.dias_peak = nan(1,L);      positions_out.dias_peak(index_R) = positions.dias_peak;

positions = positions_out;


end

