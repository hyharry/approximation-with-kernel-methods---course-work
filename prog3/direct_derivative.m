function [div_a_grad_k, a_grad_k] = direct_derivative(a, grad_k, grad_a, laplace_k)
% DIRECT_DERIVATIVE use laplace of kernel to construct function handle
% 
% a: func_handle diffusivity
% grad_k: cell array of func_handle {grad_k_x, grad_k_y, ...}
% grad_a: gradient of diffusivity
% laplace_k: laplace of kernel function
% 
% return: div_a_grad_k: func_handle div(a*grad(k))
%         a_grad_k: cell of func_handel {a*grad(k)_x, a*grad(k)_y, ...}

% $Author: yihu $	$Date: 2016/01/21 17:14:12 $	$Revision: 0.1 $

dim = length(grad_k);
a_grad_k = cell(1,dim);

grad_a_grad_k = cell(1,dim);

% Generat function handle for each axis
for i=1:dim
    % a*grad(k) in current direction
    a_grad_k{i} = @(Xcol, Xcenter) bsxfun(@times, a(Xcol), grad_k{i}(Xcol, Xcenter));
    subindex = @(A,c) A(:,c); 
    grad_a_grad_k{i} = @(Xcol, Xcenter) bsxfun(@times, subindex(grad_a(Xcol),i), grad_k{i}(Xcol, Xcenter));
end

a_laplace_k = @(Xcol, Xcenter) bsxfun(@times, a(Xcol), laplace_k(Xcol, Xcenter));
grad_a_grad_k_sum = @(Xcol, Xcenter) sum_func_handle(Xcol, Xcenter, grad_a_grad_k);

div_a_grad_k = @(Xcol, Xcenter) sum_func_handle(Xcol, Xcenter, {grad_a_grad_k_sum, a_laplace_k});
end