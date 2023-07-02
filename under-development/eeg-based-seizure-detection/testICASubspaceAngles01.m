close all;
clear;
clc;

% parameters
subject = {'Dog_1', 'Dog_2', 'Dog_3', 'Dog_4', 'Dog_5', 'Patient_1', 'Patient_2'};
ffs = 160; % resample rate
f0 = 60; % powerline frequency Hz
Q = 30; % notch filter Q-factor
epochlen = 60; % epoch length in seconds
epochoverlap = 50; % overlap between epochs in seconds
sourcenum = 6;
resultfname = 'testICASubspaceAngles01.txt';

x = [];
for m = 1:length(subject)
    path = ['E:\Sameni\Projects\Seizure\' subject{m} '\'];
    d = dir([path '*.mat']);
    NumRecords = length(d);
    for k = 1 : NumRecords,
        if(~isempty(x))
            clear x;
        end
        load([path d(k).name]);
        reg = regexp(d(k).name, '_');
        varname = [d(k).name(reg(2)+1:reg(4)-1) '_' num2str(str2double(d(k).name(end-7:end-4)))];
        md = d(k).name(reg(2)+1:reg(3)-1);
        if(isequal(md, 'interictal'))
            mode = 1;
        elseif(isequal(md, 'preictal'))
            mode = 2;
        elseif(isequal(md, 'test'))
            mode = 3;
        end
        
        eval(['data = ' varname '; clear ' varname ';']);
        
        x = data.data;
        fs = data.sampling_frequency;
        clear data;
        
        % resample data
        if(fs > ffs) % for the patients who have high sampling rates
            x = resample(x', ffs, round(fs))';
            fs = ffs;
        end
        
        % pre-process
        % x = LPFilter(x, 80.0/fs); % lowpass filter
        x = x - LPFilter(x, 0.1/fs); % highpass filter
        % IIR notch
        Wo = f0/(fs/2);  BW = Wo/Q;
        [b,a] = iirnotch(Wo,BW);
        x = filter(b, a, x, [], 2);
        
        % normalize
        x = (x - mean(x, 2)*ones(1, size(x,2)))./(std(x, [], 2)*ones(1, size(x,2)));
        % PlotECG(x, 4, 'b', fs);
        
        channelnum = size(x,1);
        epochnum = floor((length(x)-epochlen*fs)/((epochlen-epochoverlap)*fs));
        
        % BSS
        wxx = zeros(epochnum, sourcenum, channelnum);
        wPCA = zeros(epochnum, sourcenum, channelnum);
        wJADE = zeros(epochnum, sourcenum, channelnum);
        %         wSOBI = zeros(epochnum, sourcenum, channelnum);
        %         wSCA = zeros(epochnum, sourcenum, channelnum);
        for i = 1:epochnum,
            xx = x(:, round(i*(epochlen-epochoverlap)*fs:i*(epochlen-epochoverlap)*fs+epochlen*fs));
            
            % PCA
            Cx = cov(xx');
            [V D] = eig(Cx);
            [Y,I] = sort(diag(D), 1, 'descend');
            wPCA(i,:,:) = V(:, I(1:sourcenum))';
            
            % JADE
            wJADE(i,:,:) = jadeR(xx, sourcenum);
            wxx(i,:,:) = normr(squeeze(wJADE(i,:,:)));
            
            % SOBI
            %             w = pinv(sobi(xx));
            %             wSOBI(i,:,:) = w(1:sourcenum, :);
            
            % SCA
            %             [~, w] = SCA2(xx, 10.0/fs, 10.0/fs, 5);
            %             wSCA(i,:,:) = w(1:sourcenum, :);
            
            %     figure
            %     plot(xx');
            %     grid
            %     title(num2str(i));
        end
        
        % Subspace robustness and angles
        prodsx = zeros(epochnum, epochnum, sourcenum, sourcenum);
        prods2x2x = zeros(epochnum*sourcenum, epochnum*sourcenum);
        averageonepochsx = zeros(epochnum*sourcenum, 1);
        epochcoherencex = zeros(epochnum, epochnum);
        AnglesPCA = zeros(epochnum, 1);
        AnglesJADE = zeros(epochnum, 1);
        %         AnglesSOBI = zeros(epochnum, 1);
        %         AnglesSCA = zeros(epochnum, 1);
        for i = 1:epochnum, % over epochs
            AnglesPCA(i) = subspace(squeeze(wPCA(i,:,:))', squeeze(wPCA(1,:,:))')*180/pi;
            AnglesJADE(i) = subspace(squeeze(wJADE(i,:,:))', squeeze(wJADE(1,:,:))')*180/pi;
            %             AnglesSOBI(i) = subspace(squeeze(wSOBI(i,:,:))', squeeze(wSOBI(1,:,:))')*180/pi;
            %             AnglesSCA(i) = subspace(squeeze(wSCA(i,:,:))', squeeze(wSCA(1,:,:))')*180/pi;
            for j = 1:epochnum, % over epochs
                prodsx(i, j, :, :) = abs(squeeze(wxx(i,:,:)) * squeeze(wxx(j,:,:))');
                prx = squeeze(prodsx(i, j, :, :));
                epochcoherencex(i, j) = mean(mean(prx));
                prods2x2x((i-1)*sourcenum + 1 : i*sourcenum, (j-1)*sourcenum + 1 : j*sourcenum) = prx;
                if(i~=j)
                    colmaxx = max(prx, [], 2);
                    averageonepochsx((i-1)*sourcenum + 1 : i*sourcenum, :) = averageonepochsx((i-1)*sourcenum + 1 : i*sourcenum, :) + colmaxx;
                end
            end
        end
        averageonepochsx = averageonepochsx/(epochnum - 1);
        
        f1 = median(epochcoherencex(:));
        f2 = mean(epochcoherencex(:));
        f3 = std(epochcoherencex(:));
        
        f4 = median(prods2x2x(:));
        f5 = mean(prods2x2x(:));
        f6 = std(prods2x2x(:));
        
        f7 = median(averageonepochsx(:));
        f8 = mean(averageonepochsx(:));
        f9 = std(averageonepochsx(:));
        
        % average PCA subspace angles
        f10 = median(AnglesPCA(:));
        f11 = mean(AnglesPCA(:));
        f12 = std(AnglesPCA(:));
        
        % average JADE subspace angles
        f13 = median(AnglesJADE(:));
        f14 = mean(AnglesJADE(:));
        f15 = std(AnglesJADE(:));
        
        %         % average SOBI subspace angles
        %         f16 = median(AnglesSOBI(:));
        %         f17 = mean(AnglesSOBI(:));
        %         f18 = std(AnglesSOBI(:));
        %
        %         % average SCA subspace angles
        %         f19 = median(AnglesSCA(:));
        %         f20 = mean(AnglesSCA(:));
        %         f21 = std(AnglesSCA(:));
        
        % write results
        fid = fopen(resultfname,'a');
        %         fprintf(fid, '%s\t%d\t%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n', d(k).name, m, k, mode, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20, f21);
        fprintf(fid, '%s\t%d\t%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n', d(k).name, m, k, mode, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15);
        fclose(fid);
        disp(['subject: ', d(k).name]);
        %         end % end if mode(k) == 3
    end
    % % %     pack
end

clock