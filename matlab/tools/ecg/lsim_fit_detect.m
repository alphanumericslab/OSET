function viterbi_trans_detection = lsim_fit_detect(obs_seqment_beats, type_det, init_state_1, init_state_2)

init_state_1 = cell2mat(init_state_1);
init_state_2 = cell2mat(init_state_2);


C=2;
max_state = 2;

trans_temp = [10^5,1;0,1];
trans_temp = trans_temp./sum(trans_temp,2);

trans_temp1 = [1,0;0,1];
trans_temp1 = trans_temp1./sum(trans_temp1,2);

init_model.transition_matrices_IM_cell{1,1} = trans_temp;
init_model.transition_matrices_IM_cell{2,1} = trans_temp1;
init_model.transition_matrices_IM_cell{1,2} = trans_temp1;
init_model.transition_matrices_IM_cell{2,2} = trans_temp;


temp_pi0 = [1;zeros(max_state-1,1)];
init_model.pi_0_lsim{1,1} = temp_pi0/sum(temp_pi0);
init_model.pi_0_lsim{2,1} = temp_pi0/sum(temp_pi0);

init_model.coupling_tetha_IM = [0.99,0.01;0.01,0.99];

extra.sigma_diag = 1;
extra.plot = 0;
extra.check_convergence = 0;
extra.time_series = 0;
extra.left2right = 1;
max_itration = 30;
batch_size = 30;
num_gmm_component = 1*ones(1,C);
channel_num_states = max_state*ones(1,C);


gmm_para_ch1(1).P = 1;
gmm_para_ch1(1).sigma(1).x = var(init_state_1(1,:));
gmm_para_ch1(1).mu(1).x = mean(init_state_1(1,:));

gmm_para_ch2(1).P = 1;
gmm_para_ch2(1).sigma(1).x = var(init_state_1(2,:));
gmm_para_ch2(1).mu(1).x = mean(init_state_1(2,:));

gmm_para_ch1(2).P = 1;
gmm_para_ch1(2).sigma(1).x = var(init_state_2(1,:));
gmm_para_ch1(2).mu(1).x = mean(init_state_2(1,:));

gmm_para_ch2(2).P = 1;
gmm_para_ch2(2).sigma(1).x = var(init_state_2(2,:));
gmm_para_ch2(2).mu(1).x = mean(init_state_2(2,:));

gmm_para_lsim_cell{1,1}.gmm_para = gmm_para_ch1;
gmm_para_lsim_cell{2,1}.gmm_para = gmm_para_ch2;
init_model.gmm_para_lsim = gmm_para_lsim_cell;


num_batch = max(1,floor(size(obs_seqment_beats,2)/batch_size));


parfor b = 1:num_batch

    if b~=num_batch
        channels_observations = obs_seqment_beats(:,(b-1)*batch_size+1:b*batch_size);
    else
        channels_observations = obs_seqment_beats(:,(b-1)*batch_size+1:size(obs_seqment_beats,2));
    end

    [~ , coupling_tetha_convex_comb , transition_matrices_convex_comb ,  lsim_gmm_para ] =  em_lsim( channels_observations , channel_num_states , num_gmm_component , max_itration , extra, init_model);
    [pi_0_ehmm{b}, coupling_tetha_ehmm{b},  transition_ehmm{b}, ehmm_gmm_para{b}] = im_para_eqhmm(init_model.pi_0_lsim, lsim_gmm_para, coupling_tetha_convex_comb, transition_matrices_convex_comb);

end

viterbi_trans_detection = nan(size(obs_seqment_beats,2),2);

for p = 1:size(obs_seqment_beats,2)

    b = floor(p/batch_size)+1;
    b(b>num_batch) = num_batch;

    temp = obs_seqment_beats(:,p);
    Y_in{1} = cell2mat(temp);

    try
        [~ , X_star]  = viterbi_chmm( pi_0_ehmm{b} , coupling_tetha_ehmm{b} ,  transition_ehmm{b}  , ehmm_gmm_para{b},  Y_in  );

        if strcmp(type_det,'on')
            viterbi_trans_detection(p,1) = find(X_star==4,1,"first");
            viterbi_trans_detection(p,2) = find(X_star==1,1,"last")+1;
        elseif  strcmp(type_det,'off')
            viterbi_trans_detection(p,1) =  find(X_star==4,1,"first")-1;
            viterbi_trans_detection(p,2) =  find(X_star==1,1,"last");
        end

    catch
        try
            channels_observations = temp;
            [~ , coupling_tetha_convex_comb , transition_matrices_convex_comb ,  lsim_gmm_para ] =  em_lsim( channels_observations , channel_num_states , num_gmm_component , 10 , extra, init_model);
            [pi_0_ehmm_temp, coupling_tetha_ehmm_temp,  transition_ehmm_temp, ehmm_gmm_para_temp] = im_para_eqhmm(init_model.pi_0_lsim, lsim_gmm_para, coupling_tetha_convex_comb, transition_matrices_convex_comb);

            [~ , X_star]  = viterbi_chmm( pi_0_ehmm_temp, coupling_tetha_ehmm_temp,  transition_ehmm_temp, ehmm_gmm_para_temp,  Y_in  );

            if strcmp(type_det,'on')
                viterbi_trans_detection(p,1) = find(X_star==4,1,"first");
                viterbi_trans_detection(p,2) = find(X_star==1,1,"last")+1;
            elseif  strcmp(type_det,'off')
                viterbi_trans_detection(p,1) =  find(X_star==4,1,"first")-1;
                viterbi_trans_detection(p,2) =  find(X_star==1,1,"last");
            end

        catch
            viterbi_trans_detection(p,:) =   viterbi_trans_detection(p,:) ;
        end

    end


end




