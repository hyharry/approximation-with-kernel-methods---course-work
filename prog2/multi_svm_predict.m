function y_res = multi_svm_predict(Xtest,X,y,k,svm_array,cls_label)

if isrow(y)
    y = y(:);
end

svm_num = length(svm_array);
ytest = zeros(svm_num,size(Xtest,2));

for i = 1:svm_num
    % select the corresponding data to bin_svm
    [X_i, ~] = multi_to_binary(X,y,cls_label(i,1),cls_label(i,2));
    Ktest = k(Xtest, X_i);
    
    % predict for ith svm
    [ytest(i,:), ~] = svm_array{i}(Ktest);
end

% get the most occurence for each test data
y_res = mode(ytest);