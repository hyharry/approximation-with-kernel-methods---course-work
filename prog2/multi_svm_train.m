function [svm_array,cls_label] = multi_svm_train(X,y,k,C)
% svm_array = {svm_func_handle1, svm_func_handle2, ...}
% cls_label = [svm1_cls1, svm1_cls2;
%              svm2_cls1, svm2_cls2; ...]

label = unique(y);
nc = length(label);
svm_num = nc*(nc-1)/2;
svm_array = cell(1,svm_num);
cls_label = zeros(svm_num,2);

count = 1;
for i = 1:nc
    for j = i+1:nc
        % select cls1 and cls2, and transform the data for bin_svm
        [X_train, y_train] = multi_to_binary(X,y,i,j);
        K_train = k(X_train, X_train);
        [alpha, b] = binary_svm_train(K_train,y_train,C);
        
        % generate bin_svm based on cls1 and cls2
        svm_array{count} = @(Ktest) binary_svm_predict(Ktest,y_train,alpha,b,[i,j]);
        cls_label(count,:) = [i,j];
%         display(sprintf('i = %i, j = %i, training completed', i, j));
        count = count + 1;
    end  
end
