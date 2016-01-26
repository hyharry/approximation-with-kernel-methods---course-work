function [A, b] = assemble_RBF_collocation_system(model, kernel, kernel_grad, Xcol, Xcenter, ... 
                                                    fd_para, kernel_laplace)
%ASSEMBLE_RBF_COLLOCATION_SYSTEM: assemble system matrix and rhs for general
%elliptic PDE. Difference quotients are used when assembling for the derivatives
%
% [A, b] = ASSEMBLE_RBF_COLLOCATION_SYSTEM(model, kernel, kernel_grad, Xcol, Xcenter)
% 
% model: see example poisson_model.m
% Xcenter, Xcol: coord of center and col, each row for a node
% kernel: func_handle return kernel matrix, like previous prog exercise
% kernel_grad: func_handle_cell {grad_x_matrix, grad_y_matrix, ...} = my_kernel_grad
% kernel_laplace: alternative formulation with kernel laplace
% fd_para: finite difference parameter
% kernel_laplace: function to generate kernel laplace matrix

% $Author: yihu $	$Date: 2016/01/21 15:00:30 $	$Revision: 0.1 $


% Group points into inner, dirichlet, neumman
node_type = model.boundary_type(Xcol);
inner_ind = node_type == 0;
X_inner = Xcol(inner_ind,:);
diri_ind = node_type == -1;
X_diri = Xcol(diri_ind,:);
neum_ind = node_type == -2;
X_neum = Xcol(neum_ind,:);

% Equation (1), PDE in domain
a = model.diffusivity;

display('func_handle time')
tic;
% Generate function handle for the discretized equation
% TODO: routine that calculates result matrices is more efficient
if nargin == 7
    a_grad = model.diffusivity_gradient;
    [div_a_grad_k, a_grad_k] = direct_derivative(a, kernel_grad, a_grad, kernel_laplace);
else
    [div_a_grad_k, a_grad_k] = diff_quotients(a, kernel_grad, fd_para);
end
toc;

display('assemble time')
tic;
A1 = -div_a_grad_k(X_inner, Xcenter);
b1 = model.source(X_inner);
toc;

% Equation (2), Dirichlet BC
A2 = kernel(X_diri, Xcenter);
b2 = model.dirichlet_values(X_diri);

% Equation(3), Neumann BC
normal = model.normals(X_neum);
A3 = bsxfun(@times, normal(:,1), a_grad_k{1}(X_neum, Xcenter)) + ...
    bsxfun(@times, normal(:,2), a_grad_k{2}(X_neum, Xcenter));
b3 = model.neumann_values(X_neum);

A = [A1;A2;A3];
b = [b1;b2;b3];

% A = [A1;A2];
% b = [b1;b2];

end
