function classifier(Xtest, X, y, C, k)
% binary classifier 

K_train = k(X,X);
K_test = k(Xtest,X);
[alpha, b] = binary_svm_train(K_train,y,C);
[ytest, ftest] = binary_svm_predict(K_test, y, alpha, b);

% visualization
n = sqrt(size(Xtest,2));
Xmesh = reshape(Xtest(1,:), [n,n]);
Ymesh = reshape(Xtest(2,:), [n,n]);
Zmesh = reshape(ftest,[n,n]);

% [C,h] = contour(Xmesh,Ymesh,Zmesh);
[C,h] = contourf(Xmesh,Ymesh,Zmesh);
h.LevelStep = 1;
% clabel(C,h);
v = [-1,0,1];
clabel(C,h,v);

hold on
scatter(X(1,:), X(2,:),'r','filled')
% gscatter(Xtest(1,:)', Xtest(2,:)',ytest');
