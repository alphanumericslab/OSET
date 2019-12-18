function score = FSCost(XT, yT, Xt, yt)

svmStruct = svmtrain(XT, yT, 'kernel_function', 'rbf', 'rbf_sigma', 100, 'boxconstraint', 5);
% svmStruct = svmtrain(r_train, types(train_index_subset), 'kernel_function', 'mlp');%, 'rbf_sigma', rbf_sigma, 'boxconstraint', length(interictal_train_index_subset)/length(preictal_train_index_subset));
estimated_class = svmclassify(svmStruct, Xt);
score = sum(estimated_class == yt);

end