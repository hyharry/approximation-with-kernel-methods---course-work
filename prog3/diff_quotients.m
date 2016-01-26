function [div_a_grad_k, a_grad_k] = diff_quotients(a, grad_k, para)
% DIFF_QUOTINETS generate function handle for later assembling using pure
% functional programming. Finite quotient of derivatice approximation is used
% 
% a: func_handle diffusivity
% grad_k: cell array of func_handle {grad_k_x, grad_k_y, ...}
% para: fd distance
% 
% return: div_a_grad_k: func_handle div(a*grad(k))
%         a_grad_k: cell of func_handel {a*grad(k)_x, a*grad(k)_y, ...}

% $Author: yihu $	$Date: 2016/01/21 17:14:12 $	$Revision: 0.1 $

dim = length(grad_k);
a_grad_k = cell(1,dim);


diff_a_grad_k = cell(1,dim);
h = zeros(1, dim);

% Generat function handle for each axis
for i=1:dim
    % Difference for current direction
    h(i) = para;

    % Function handle to generate X_plus and X_minus
    X_plus = @(X_inner) bsxfun(@plus, X_inner, h);
    X_minus = @(X_inner) bsxfun(@minus, X_inner, h);

    % a*grad(k) in current direction
    a_grad_k{i} = @(Xcol, Xcenter) bsxfun(@times, a(Xcol), grad_k{i}(Xcol, Xcenter));

    % Finite Difference approximation of divergence
    a_grad_k_plus = @(Xcol, Xcenter) a_grad_k{i}(X_plus(Xcol), Xcenter);
    a_grad_k_minus = @(Xcol, Xcenter) a_grad_k{i}(X_minus(Xcol), Xcenter);
    diff_a_grad_k{i} = @(Xcol, Xcenter) (a_grad_k_plus(Xcol, Xcenter) - a_grad_k_minus(Xcol, Xcenter))/(2*para);
end

div_a_grad_k = @(Xcol, Xcenter) sum_func_handle(Xcol, Xcenter, diff_a_grad_k);

end

