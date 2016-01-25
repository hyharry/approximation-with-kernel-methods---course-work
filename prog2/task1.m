function task1

% test set
Xtest = test_point_gen(-1,2,40);

% c) Application on linearly separable data
load('data2d_linear.mat','-mat');
C = [0.2,2,10,100];
n = length(C);
k = @(X1,X2) k_linear(X1,X2);
gcf
for i=1:n
    subplot(2,2,i);
    classifier(Xtest,X,y,C(i),k);
    title(sprintf('C = %f',C(i)));
end

% d) Application on nonlinearly separable data
load('data2d_nonlinear.mat','-mat');
gamma = [0.2,2,10,100];
n = length(gamma);
figure
for i=1:n
    k = @(X1,X2) k_gauss(X1,X2,gamma(i));
    subplot(2,2,i);
    classifier(Xtest,X,y,Inf,k);
    title(['\gamma =', num2str(gamma(i))]);
end