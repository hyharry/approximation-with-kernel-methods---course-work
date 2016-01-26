function [error_l] = Task3(alpha, m, gamma_list, method)
%TASK3 Summary of this function goes here
% 
% [error] = TASK3(m, gamma_list, method) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: yihu $	$Date: 2016/01/21 23:04:06 $	$Revision: 0.1 $

refine = 2; % for visualization

beta = 0.2;
model = general_model([alpha beta]);

a = -1;
b = 1;
x = linspace(a,b,m);
[X1, X2] = meshgrid(x,x);
Xcenter = [X1(:) X2(:)];
% [Xcenter, q] = node_sample(100, 10,-1,1,-1,1);
row_ind = Xcenter(:,1)>0 & Xcenter(:,2)>0;
% bound_ind_1 = Xcenter(:,1)<-1+1e-10 & Xcenter(:,2)<=-1+beta;
% bound_ind_2 = Xcenter(:,2)<-1+1e-10 & Xcenter(:,1)<=-1+beta;
% row_ind = row_ind | bound_ind_2 | bound_ind_1;
Xcenter(row_ind,:) = [];

xvis = linspace(a,b,refine*m);
[Xvis, Yvis] = meshgrid(xvis,xvis);
ind = Xvis>0 & Yvis>0;
Xvis(ind) = NaN;
Yvis(ind) = NaN;

% Exclude upper right corner
% row_ind = Xvis(:,1)>0 & Xvis(:,2)>0;
% Xvis(row_ind,:) = [];

% corners = [-1, -1;
%            -1, 1;
%            0, 1;
%            0, 0;
%            1, 0;
%            1, -1];
% Xcenter = trim_data_set(Xcenter, corners);
% Xvis = trim_data_set(Xvis, corners);  

test_num = length(gamma_list);
f_approx = cell(1,test_num);

if method == 'd'
    % Assemble with finite quotient
    for i=1:length(gamma_list)
        gamma = gamma_list(i);
        kernel = @(X1,X2) k_gauss(X1,X2,gamma);
        kernel_grad = k_gauss_grad(2,gamma);
        fd_para = 0.01;

        assemble_method = @(model,kernel,Xcol,Xcenter) ... 
                            assemble_RBF_collocation_system(model, kernel, kernel_grad, ... 
                                                            Xcol, Xcenter, fd_para);

        f_approx{i} = solve_pde(Xcenter, model, kernel, assemble_method);
    end
elseif method == 'f'
    % Assemble with laplace of kernel
    for i=1:length(gamma_list)
        gamma = gamma_list(i);
        kernel = @(X1,X2) k_gauss(X1,X2,gamma);
        kernel_grad = k_gauss_grad(2,gamma);
        fd_para = 0.01;

        kernel_laplace = k_gauss_laplace(gamma);
        assemble_method = @(model,kernel,Xcol,Xcenter) ... 
                            assemble_RBF_collocation_system(model, kernel, kernel_grad, Xcol, Xcenter, ... 
                                                            fd_para, kernel_laplace);

        f_approx{i} = solve_pde(Xcenter, model, kernel, assemble_method);
    end
else
    error('method not valid');
end

error_l = post_process(model.solution, f_approx, gamma_list, Xvis, Yvis, Xcenter);


function X = trim_data_set(X, X_exclude)
% Assistance function to exclude points in point set
for i=1:size(X_exclude,1)
    X_tmp = bsxfun(@(x,y) abs(x-y), X, X_exclude(i,:));
    ind_matrix = X_tmp<=1e-10;
    ind = prod(ind_matrix,2)==1;
    X(ind, :) = [];
end