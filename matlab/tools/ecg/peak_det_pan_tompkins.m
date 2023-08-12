function [peaks, peak_indexes] = peak_det_pan_tompkins(ecg_data, fs, onoff_line)
% peak_det_pan_tompkins - R-peak detector based on Pan-Tompkins method.
% QRS Detection with Pan-Tompkins algorithm.
% Pan et al. 1985: A Real-Time QRS Detection Algorithm.
%   [peaks, peak_indexes] = peak_det_pan_tompkins(data, fs, varargin)
%
%   This function implements the Pan-Tompkins algorithm for R-peak detection
%   in ECG signals, with a simplified post-detection R-peak selection logic
%
%   Inputs:
%       data: Vector of input ECG data
%       fs: Sampling rate in Hz
%       fc_low (Optional): BP filter lower cutoff frequency in Hz (default: 5.0 Hz)
%       fc_high (Optional): BP filter upper cutoff frequency in Hz (default: 15.0 Hz)
%       window_length (Optional): Integration window length in seconds (default: 0.150 s)
%       threshold_ratio (Optional): Threshold ratio for peak detection (default: 0.2)
%       refractory_period (Optional): Refractory period in seconds (default: 0.2 s)
%
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
%
%   Reference:
%       Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
%       Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532
%
%   Reza Sameni, 2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% (optional input):
%        - refrac:      refractory time [ms] (default 200ms)
%        - refracT:     refractory time for T-wave [ms] (default 360ms).
%                       Part of the heuristics behind the algorithm.
%        _ QRSmaxwidth: expected maximal length of QRS complexes [ms].
%                       Default 150ms.

QRSmaxwidth = 0.150;
refrac = 0.200;
refracT = 0.360;

% allow one of both thresholds I or F to be deceeded if associated region is detected   
regionFLAG = 1;
verbose = 1;
proc_pre_QRS = 0;

data = ecg_data(:);

fs_pt = 200;
total_delay = 0;

G = gcd( fs_pt,fs );
data = resample(data,fs_pt/G,fs/G);

%% Pan-Tomkins signal processing

% #1 remove the mean
data = data - mean(data);

% #2 low-pass filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a = [1 -2 1];
filtered_data = filter(b,a,data);
total_delay = total_delay + 6;

% #3 high-pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
b = zeros(1,33);
b(1) = -1; b(17) = 32; b(33) = 1;
a = [1 1];
filtered_data = filter(b,a,filtered_data);    % Without Delay
total_delay = total_delay + 16;
filter_delay = total_delay;

% Differentiation to enhance R-peaks
% H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
b = [1 2 0 -2 -1].*(1/8)*fs;

diff_data = filter(b,1,filtered_data);
total_delay = total_delay + 2;

% Squaring to further emphasize R-peaks
squared_data = diff_data .^ 2;

% Moving average integration
window_length = round( QRSmaxwidth * fs_pt/2);
integrated_data = movmean(squared_data, [window_length,window_length]);

% remove delays & resample
diff_data = [diff_data(total_delay:end); repelem(diff_data(end),total_delay-1,1)]';
diff_data = resample(diff_data,fs/G ,fs_pt/G);

integrated_data = [integrated_data(total_delay:end); repelem(integrated_data(end),total_delay-1,1)]';
integrated_data = resample(integrated_data,fs/G ,fs_pt/G);


filtered_data = [filtered_data(filter_delay:end); repelem(filtered_data(end),filter_delay-1,1)]';
filtered_data = resample(filtered_data,fs/G ,fs_pt/G);

data = ecg_data(:)';

%% Adaptive thresholding

refrac = round( refrac * fs);
refracT = round( refracT * fs);

%=======================================================================
% ===================  learning phase 1 ==================================
%=======================================================================


% initializing thresholds based on first 2s of signal
peakI = max(integrated_data(1:(2*fs)));
thresholdI1 = 0.125*peakI;
%thresholdI2 = 0.5*thresholdI1; % is automatically updated if necessary
peakF = max(filtered_data(1:(2*fs)));
thresholdF1 = 0.125*peakF;
%thresholdF2 = 0.5*thresholdF1; % is automatically updated if necessary

%=======================================================================
% ========================= learning phase 2 ==============================
%=======================================================================

% (this isn't explicitly specified in the originally proposed paper)
learned = 0;
RR_init = 0;
i = 1;

%including possibility to reduce threshold to one half for five times during learning phase 2
%determine threshold-exceeding signal excerpts
while ((learned == 0)&&(i <= 6))
    flagS = 0;
    flagE = 0;

    aboveI1 = (integrated_data>thresholdI1);
    aboveF1 = (filtered_data>thresholdF1);
    %exceeding both thresholds
    aboveCommon = aboveI1.*aboveF1;
    negativeCommon = ~(aboveCommon);
    %assign regions
    aboveCommonRegion = aboveCommon;
    found_region = 0;
    start_index = 0;
    for index = 1:length(aboveCommonRegion)
        %starting point of region
        if(aboveCommonRegion(1,index) == 1)
            found_region = 1;
            start_index = index;
        end
        %assign region
        if(found_region)
            aboveCommonRegion(1,index) = aboveI1(1,index)|aboveF1(1,index);
            %detect end of region
            regI1 = diff(aboveI1(1,start_index:index));
            regF1 = diff(aboveF1(1,start_index:index));
            if(~((sum(abs(regI1))==0)||(sum(abs(regF1))==0)))
                found_region = 0;
            end
            if((sum(abs(regI1))>0)&&(sum(abs(regF1))>0))
                aboveCommonRegion(1,start_index+1:index) = 0;
            end
        end
    end

    %count exceeding regions
    aboveBackup = aboveCommon;
    if (regionFLAG)
        aboveCommon = aboveCommonRegion;
    end

    count_thresAbove = sum(abs(diff(aboveCommon)));
    if (aboveCommon(1,1) == 1)
        count_thresAbove = count_thresAbove + 1;
        flagS = 1;
    end
    if (aboveCommon(1,length(aboveCommon)) == 1)
        count_thresAbove = count_thresAbove + 1;
        flagE = 1;
    end
    count_thresAbove = count_thresAbove/2;
    %threshold re-adjustment for possible learning phase recall (learned = 0)
    thresholdI1 = thresholdI1/2;
    thresholdF1 = thresholdF1/2;
    i = i + 1;
    if (i >= 3)&& verbose
        %disp(['QRS-detector notification (Pan&Tompkins):']);
        disp('threshold reduced during learning phase' );
    end

    % generation of corresponding index-pairs for threshold-exceeding regions
    start_region_index = (find(diff(aboveCommon) == 1))+1;
    end_region_index = find(diff(aboveCommon) == -1);
    if (flagS)
        start_region_index = [1,start_region_index];
    end
    if (flagE)
        end_region_index = [end_region_index, length(aboveCommon)];
    end

    j = 1;
    k = 2;

    %determine first two reliable initial qrs-complexes
    while ((learned == 0)&&(j < count_thresAbove)&&(k <= count_thresAbove)&&(count_thresAbove >= 2))

        %pick out current threshold-exceeding signal regions
        first_region = data(start_region_index(j):end_region_index(j));
        second_region = data(start_region_index(k):end_region_index(k));

        %calculating distance of peaks using input signal
        [cmax1,imax1] = max(first_region);
        imax1 = imax1+start_region_index(j)-1;
        [cmax2,imax2] = max(second_region);
        imax2 = imax2+start_region_index(k)-1;
        distance = imax2-imax1;

        %check refractory/t-wave criterion and manipulate indices
        %inside refractory period
        if (distance <= refrac)
            k = k+1;
            %between refractory period and t-wave-citerion
        elseif ((distance <= refracT)&&(distance > refrac))
            slope_ratio = max(diff_data(start_region_index(j):end_region_index(j)))/max(diff_data(start_region_index(k):end_region_index(k)));
            if (slope_ratio > 2)
                %probably t-wave
                k = k+1;
            elseif (slope_ratio < 0.5)
                %j = j+1; %preceeding QRS-detection is identified as t-wave
                %k = k+1;
                learned = 1;
            else
                learned = 1;
            end

        else
            learned = 1;
        end

    end

end

if (learned ~= 0)
    qrs_pos(1,1) = 0;
    qrs_pos(1,2) = imax1;
    qrs_pos(1,3) = imax2;
    %re-initialize thresholds for qrs-detection regarding first two detections
    RR_init = imax2-imax1;
    if verbose;disp(['RR_init: ' num2str(RR_init/fs) ' s']);end

    posMask = aboveBackup(start_region_index(j):end_region_index(k));
    negMask = negativeCommon(start_region_index(j):end_region_index(k));

    peakIs = max(posMask.*integrated_data(start_region_index(j):end_region_index(k)));
    peakIn = max(negMask.*integrated_data(start_region_index(j):end_region_index(k)));

    peakFs = max(posMask.*filtered_data(start_region_index(j):end_region_index(k)));
    peakFn = max(negMask.*filtered_data(start_region_index(j):end_region_index(k)));

    SPKI = peakIs;
    NPKI = peakIn;
    SPKF = peakFs;
    NPKF = peakFn;

    thresholdI1m(1,1) = thresholdI1*2;
    thresholdF1m(1,1) = thresholdF1*2;

    % [cmax, starting_detection] = max(posMask.*data(start_region_index(j):end_region_index(k)));
    % qrs_pos(1,1) = starting_detection+start_region_index(j)-1;
    qrs_pos(1,1) = imax1;
end


%=======================================================================
% =====================  detection phase ===============================
%=======================================================================

%initialization of detection phase
qrs_pos = qrs_pos(1,1);

RR_AV1vec = NaN(1,8);
RR_AV1vec(1,1) = RR_init;
RR1 = 0;

RR_AV2vec = NaN(1,8);
RR_AV2vec(1,1) = RR_init;

finished = 0;
thres_add = 1;
end_add = 0;
segment_OK = 0;
initial_TA = 1;
proc_per_QRS = 1;
actual_index = qrs_pos(1,1);

while ((finished == 0)&&(ceil((1.66*RR1)+actual_index)<=length(data)))

    if(segment_OK)
        end_add = 0;
        thres_add = 1;
        initial_TA = 1;
    end

    RR1 = mean(RR_AV1vec(find(~isnan(RR_AV1vec))));
    RR2 = mean(RR_AV2vec(find(~isnan(RR_AV2vec))));

    %regulary heart rate check
    if (RR1 == RR2)
        factor = 1;
    else
        factor = 0.5;
    end

    %updating thresholds
    thresholdI1 = factor*(NPKI + (0.25*(SPKI-NPKI)))*thres_add;
    %thresholdI2 = 0.5*thresholdI1; %is automatically updated if necessary

    thresholdF1 = factor*(NPKF + (0.25*(SPKF-NPKF)))*thres_add;
    %thresholdF2 = 0.5*thresholdF1; %is automatically updated if necessary

    %begin new search from actual detection
    %building up data and apply threshold
    if ((actual_index+(1.66*RR1)+end_add) <= length(data))
        end_segment = ceil((actual_index+(1.66*RR1)))+end_add;
        actual_mov_dat = integrated_data(actual_index:end_segment);
        actual_filt_dat = filtered_data(actual_index:end_segment);
        actual_diff_dat = diff_data(actual_index:end_segment);

    else
        break;
    end

    i = 1;
    segment_OK = 0;


    %including possibility to reduce threshold to one half for one time
    %during detection phase of one single segment
    %determine threshold-exceeding signal excerpts
    while((segment_OK == 0)&&(i <= 2))

        flagS = 0;
        flagE = 0;

        aboveI1 = (actual_mov_dat(1,:)>thresholdI1);
        aboveF1 = (actual_filt_dat(1,:)>thresholdF1);
        %special case of every feature signal value is below threshold after threshold-halving
        if (min(actual_mov_dat(1,:)) > thresholdI1)

            actual_range = max(actual_mov_dat(1,:)) - min(actual_mov_dat(1,:));

            new_thresM = exp(-(proc_per_QRS/3));
            new_thresM = new_thresM*actual_range;
            new_thresM = min(actual_mov_dat(1,:))+new_thresM;

            aboveI1 = (actual_mov_dat(1,:) > new_thresM);

        end
        if (min(actual_filt_dat(1,:)) > thresholdF1)

            actual_range = max(actual_filt_dat(1,:)) - min(actual_filt_dat(1,:));

            new_thresF = exp(-(proc_pre_QRS/3));
            new_thresF = new_thresF*actual_range;
            new_thresF = min(actual_filt_dat(1,:))+new_thresF;

            aboveF1 = (actual_filt_dat(1,:) > new_thresF);

        end


        %exceeding both thresholds
        aboveCommon = aboveI1.*aboveF1;
        negativeCommon = ~(aboveCommon);
        %assign regions
        aboveCommonRegion = aboveCommon;
        found_region = 0;
        start_index = 0;
        for index = 1:length(aboveCommonRegion)
            %starting point of region
            if(aboveCommonRegion(1,index) == 1)
                found_region = 1;
                start_index = index;
            end
            %assign region
            if(found_region)
                aboveCommonRegion(1,index) = aboveI1(1,index)|aboveF1(1,index);
                %detect end of region
                regI1 = diff(aboveI1(1,start_index:index));
                regF1 = diff(aboveF1(1,start_index:index));
                if(~((sum(abs(regI1))==0)||(sum(abs(regF1))==0)))
                    found_region = 0;
                end
                if((sum(abs(regI1))>0)&&(sum(abs(regF1))>0))
                    aboveCommonRegion(1,start_index+1:index) = 0;
                end
            end
        end

        %treating case of beeing unable to detect most recent detection again after updating threshold
        if (aboveCommon(1,1) == 0)
            aboveCommon(1,1) = 1;
            aboveCommonRegion(1,1) = 1;
        end


        %count exceeding regions
        aboveBackup = aboveCommon;
        if (regionFLAG)
            aboveCommon = aboveCommonRegion;
        end

        count_thresAbove = sum(abs(diff(aboveCommon)));
        if (aboveCommon(1,1) == 1)
            count_thresAbove = count_thresAbove + 1;
            flagS = 1;
        end
        if (aboveCommon(1,length(aboveCommon)) == 1)
            count_thresAbove = count_thresAbove + 1;
            flagE = 1;
        end
        count_thresAbove = count_thresAbove/2;
        %threshold re-adjustment for possible searchback (segment_OK = 0)
        thresholdI1 = thresholdI1/2;
        thresholdF1 = thresholdF1/2;
        i = i + 1;
        if (i >= 3)
            proc_per_QRS = proc_per_QRS+1;
            if (verbose)
                disp(['threshold reduced during detection phase at segment ' num2str(actual_index) '-' num2str(end_segment)]);
            end
        end

        % generation of corresponding index-pairs for threshold-exceeding regions
        start_region_index = (find(diff(aboveCommon) == 1))+1;
        end_region_index = find(diff(aboveCommon) == -1);
        if (flagS)
            start_region_index = [1,start_region_index];
        end
        if (flagE)
            end_region_index = [end_region_index, length(aboveCommon)];
        end

        j = 1;
        k = 2;
        replaceQRS = 0;
        %determine the next reasonable QRS complex
        while ((segment_OK == 0)&&(k <= count_thresAbove)&&(count_thresAbove >= 2))

            %pick out current threshold-exceeding signal regions
            first_region = data((start_region_index(j)+actual_index-1):(end_region_index(j)+actual_index-1));
            second_region = data((start_region_index(k)+actual_index-1):(end_region_index(k)+actual_index-1));

            %calculating distance of peaks using input signal
            [cmax1,imax1] = max(first_region);
            imax1 = imax1+start_region_index(j)-1;
            [cmax2,imax2] = max(second_region);
            imax2 = imax2+start_region_index(k)-1;
            distance = imax2-imax1;

            %check refractory/t-wave criterion and manipulate indices
            %inside refractory period
            if (distance <= refrac)

                aboveCommon(start_region_index(k):end_region_index(k)) = 0; %due to subsequent threshold adjustment
                negativeCommon(start_region_index(k):end_region_index(k)) = 1;

                k = k+1;

                %between refractory period and t-wave-criterion
            elseif ((distance <= refracT)&&(distance > refrac))
                slope_ratio = max(actual_diff_dat(start_region_index(j):end_region_index(j)))/max(actual_diff_dat(start_region_index(k):end_region_index(k)));
                if (slope_ratio > 2)
                    %probably t-wave

                    aboveCommon(start_region_index(k):end_region_index(k)) = 0; %due to subsequent threshold adjustment
                    negativeCommon(start_region_index(k):end_region_index(k)) = 1;

                    k = k+1;

                elseif (slope_ratio < 0.5)
                    %replaceQRS = 1; %preceeding QRS-detection is identified as t-wave and replaced

                    segment_OK = 1;
                    end_add = 0;
                    thres_add = 1;
                    proc_per_QRS = 1;
                else
                    segment_OK = 1;
                    end_add = 0;
                    thres_add = 1;
                    proc_per_QRS = 1;
                end

            else
                segment_OK = 1;
                end_add = 0;
                thres_add = 1;
                proc_per_QRS = 1;

                %check if preceeding detection could be a t-wave
                if (length(qrs_pos) >= 2)
                    rr_new = (imax2+actual_index-1)-qrs_pos(1,length(qrs_pos));
                    rr_old = qrs_pos(1,length(qrs_pos))-qrs_pos(1,length(qrs_pos)-1);
                    comb_detection = rr_new+rr_old;
                    if ((comb_detection >=  0.92*RR2)&&(comb_detection <= 1.16*RR2))
                        replaceQRS = 1;
                    end
                end

            end

        end

    end

    if(segment_OK)
        %estimate detection inside input data
        detection = imax2+actual_index-1;

        %sort in new detection
        if (replaceQRS == 0)
            qrs_pos(1,1+length(qrs_pos)) = detection;

            thresholdI1m(1,1+length(thresholdI1m)) = thresholdI1*2;
            thresholdF1m(1,1+length(thresholdF1m)) = thresholdF1*2;
        else
            qrs_pos(1,length(qrs_pos)) = detection;

            thresholdI1m(1,length(thresholdI1m)) = thresholdI1*2;
            thresholdF1m(1,length(thresholdF1m)) = thresholdF1*2;
        end


        %calculate and sort in new RR-interval
        if (length(qrs_pos) >= 2)
            RR_new = qrs_pos(1,length(qrs_pos))-qrs_pos(1,length(qrs_pos)-1);

            if (replaceQRS == 0)&& verbose
                disp(['recent RR-interval: ' num2str(RR_new/fs) ' s']);
            elseif verbose
                disp(['recent RR-interval: ' num2str(RR_new/fs) ' s (replaced)']);
            end

            if (isempty(find(isnan(RR_AV1vec(1,:)))))
                if(replaceQRS == 0)
                    RR_AV1vec = RR_AV1vec(1,2:length(RR_AV1vec));
                    RR_AV1vec(1,1+length(RR_AV1vec)) = RR_new;
                else
                    RR_AV1vec(1,length(RR_AV1vec)) = RR_new;
                end
            else
                index_insert = find(isnan(RR_AV1vec(1,:)));
                if(replaceQRS == 0)
                    index_insert = index_insert(1,1);
                    RR_AV1vec(1,index_insert) = RR_new;
                else
                    index_insert = index_insert(1,1)-1;
                    RR_AV1vec(1,index_insert) = RR_new;
                end
            end

            if ((RR_new >= (0.92*RR1))&&(RR_new <= (1.16*RR1)))

                if (isempty(find(isnan(RR_AV2vec(1,:)))))
                    if(replaceQRS == 0)
                        RR_AV2vec = RR_AV2vec(1,2:length(RR_AV2vec));
                        RR_AV2vec(1,1+length(RR_AV2vec)) = RR_new;
                    else
                        RR_AV2vec(1,length(RR_AV2vec)) = RR_new;
                    end
                else
                    index_insert = find(isnan(RR_AV2vec(1,:)));
                    if(replaceQRS == 0)
                        index_insert = index_insert(1,1);
                        RR_AV2vec(1,index_insert) = RR_new;
                    else
                        index_insert = index_insert(1,1)-1;
                        RR_AV2vec(1,index_insert) = RR_new;
                    end
                end
            elseif verbose
                disp('RR-regularity dumped - out of 92-116% RR-variance');
            end
        end

        %update actual start/end index for next treated segment
        actual_index = detection;

        %update running estimates of peak-values regarding threshold usage
        if (replaceQRS == 0)
            posMask = aboveBackup(start_region_index(j):end_region_index(k));
        else
            posMask = aboveBackup(start_region_index(k):end_region_index(k));
        end
        negMask = negativeCommon(start_region_index(j):end_region_index(k));

        if (replaceQRS == 0)
            peakIs = max(posMask.*actual_mov_dat(start_region_index(j):end_region_index(k)));
            peakFs = max(posMask.*actual_filt_dat(start_region_index(j):end_region_index(k)));
        else
            peakIs = max(posMask.*actual_mov_dat(start_region_index(k):end_region_index(k)));
            peakFs = max(posMask.*actual_filt_dat(start_region_index(k):end_region_index(k)));
        end
        peakIn = max(negMask.*actual_mov_dat(start_region_index(j):end_region_index(k)));
        peakFn = max(negMask.*actual_filt_dat(start_region_index(j):end_region_index(k)));



        if (i <= 2)
            SPKI = (0.125*peakIs)+(0.875*SPKI);
            NPKI = (0.125*peakIn)+(0.875*NPKI);
            SPKF = (0.125*peakFs)+(0.875*SPKF);
            NPKF = (0.125*peakFn)+(0.875*NPKF);
        else %double speed threshold adjustment after searchback
            SPKI = (0.25*peakIs)+(0.75*SPKI);
            NPKI = (0.25*peakIn)+(0.75*NPKI);
            SPKF = (0.25*peakFs)+(0.75*SPKF);
            NPKF = (0.25*peakFn)+(0.75*NPKF);
        end

    end

    %handling case of no detection after applying both thresholds (this isn't explicitly specified in the originally proposed paper)

    if (segment_OK == 0)
        proc_per_QRS = proc_per_QRS + 1;
        if ((end_segment + end_add + fs) <= length(data))
            %adding 1s of remaining signal to the already investigated one
            end_add = end_add + fs;
            %reducing thresholds again
            thres_add = thres_add*0.5;
            if verbose;disp(['additive calculation loop due to missing QRS-detection after regular processing phase...']);end
            %compensation of threshold recalculation
            if (initial_TA)
                thres_add = thres_add*0.5;
                initial_TA = 0;
            end
        else
            finished = 1;
        end
    end
end

peaks = 0 * ecg_data;
peaks(qrs_pos) = 1;
peak_indexes = qrs_pos;


% filtered_data;
% integrated_data;
% thresholdF1m;
% thresholdI1m;
