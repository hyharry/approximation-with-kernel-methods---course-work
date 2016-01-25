function task2

load('characters.mat','-mat');

% c) Visualization 
num_per_cls = 10; % number of demos for each cls
cls_label = unique(y);
cls_num = length(cls_label);
gcf;
count = 1;
for i = 1:cls_num
    label = cls_label(i);
    ind = find(y==label, num_per_cls);
    for j = 1:num_per_cls
        subplot(10,num_per_cls,count);
        hold on
        I = mat2gray(reshape(X(:,ind(j)),[16,16]));
        imshow(I);
        count = count + 1;
    end
end

% d) Character recognition with SVM
% error = error + 1 when the classification fails for one test point
% 1. linear kernel
C_lin = [2, 20];
lin_trial_num = length(C_lin);
error_lin = zeros(1,lin_trial_num);
error_lin_train = zeros(1,lin_trial_num);
for i = 1:lin_trial_num
    k = @(X1,X2) k_linear(X1,X2);
    [svm_array,cls_label] = multi_svm_train(X,y,k,C_lin(i));
    y_res = multi_svm_predict(Xtest,X,y,k,svm_array,cls_label);
    error_lin(i) = sum(sign(abs(y_res-ytest)));
    y_res_train = multi_svm_predict(X,X,y,k,svm_array,cls_label);
    error_lin_train(i) = sum(sign(abs(y_res_train-y)));
end

% 2. gaussian kernel
C = Inf;
gamma = [0.001,0.2,2];
gauss_trial_num = length(gamma);
error_gauss = zeros(1,gauss_trial_num);
error_gauss_train = zeros(1,gauss_trial_num);
for i = 1:gauss_trial_num
    k = @(X1,X2) k_gauss(X1,X2,gamma(i));
    [svm_array,cls_label] = multi_svm_train(X,y,k,C);
    y_res = multi_svm_predict(Xtest,X,y,k,svm_array,cls_label);
    error_gauss(i) = sum(sign(abs(y_res-ytest)));
    y_res_train = multi_svm_predict(X,X,y,k,svm_array,cls_label);
    error_gauss_train(i) = sum(sign(abs(y_res_train-y)));
end

% print errors
display(sprintf(repmat('-',1,70)));
for i = 1:lin_trial_num
    display(sprintf('linear_kernel, C = %.3f, train_err = %i, test_err = %i', ...
        C_lin(i), error_lin_train(i), error_lin(i)));
end

display(sprintf(repmat('-',1,70)));
for i = 1:gauss_trial_num
    display(sprintf('gaussian_kernel, gamma = %.3f, train_err = %i, test_err = %i', ...
        gamma(i), error_gauss_train(i), error_gauss(i)));
end