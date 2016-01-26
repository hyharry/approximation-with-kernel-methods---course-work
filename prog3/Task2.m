function [error_l] = Task2(m, gamma_list, method)
%TASK2
% 
% [erro] = TASK2(m, gamma_list) 
% m: gitter resolution
% gamma_list: different gamma to test convergence

% $Author: yihu $	$Date: 2016/01/21 21:56:44 $	$Revision: 0.1 $

refine = 5; % Smoothness of surf plot
% plot_prec = 0.05; % More elements for surf plot

% Kernel center points 
x = linspace(0,1,m);
[X1,X2] = meshgrid(x,x);
Xcenter = [X1(:) X2(:)];
% [X, q] = node_sample(100, 10);

% Points for plot
x = linspace(0,1,refine*m);
[Xvis,Yvis] = meshgrid(x,x);

model = poisson_model;

test_num = length(gamma_list);
f_approx = cell(1,test_num);

if method == 'f'
    % Assemble with finite quotient
    for i=1:test_num
        gamma = gamma_list(i);
        kernel = @(X1,X2) k_gauss(X1,X2,gamma);
        kernel_grad = k_gauss_grad(2,gamma);
        fd_para = 0.01;

        % Generate assemble method using finite quotient
        assemble_method = @(model,kernel,Xcol,Xcenter) ... 
                            assemble_RBF_collocation_system(model, kernel, kernel_grad, ... 
                                                            Xcol, Xcenter, fd_para);
    
        f_approx{i} = solve_pde(Xcenter, model, kernel, assemble_method);
    end
elseif method == 'd'
    % Assemble with laplace of kernel
    for i=1:length(gamma_list)
        gamma = gamma_list(i);
        kernel = @(X1,X2) k_gauss(X1,X2,gamma);
        kernel_grad = k_gauss_grad(2,gamma);
        fd_para = 0.01;

        % Generate assemble method of direct laplace 
        kernel_laplace = k_gauss_laplace(gamma);
        assemble_method = @(model,kernel,Xcol,Xcenter) ... 
                            assemble_RBF_collocation_system(model, kernel, kernel_grad, Xcol, Xcenter, ... 
                                                            fd_para, kernel_laplace);

        f_approx{i} = solve_pde(Xcenter, model, kernel, assemble_method);
    end
else
    error('method not valid')
end

error_l = post_process(model.solution, f_approx, gamma_list, Xvis, Yvis, Xcenter);
