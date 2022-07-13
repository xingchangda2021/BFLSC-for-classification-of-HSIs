function [W]=BFLSC(data, params)
% BFLSC Training for MxN local regions

% load feret/train.mat
% load feret/train_data3.mat
% load feret/train_label3.mat
% load feret/Tra_2.mat
% data=Xtain_data; label=Train_label;

LSM =lsm_extract(data,params.M,params.N);

for i=1:params.M
    for j=1:params.N
        fprintf('learn BFLSC for LSM region (%i,%i)',i,j);
        W{i,j} = BFLSC_learn(LSM{i,j},params.w, ...
            params.n_iter,params.lambda1,params.beta,params.gamma);
    end
end

end
