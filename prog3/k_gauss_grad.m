function kernel_grad_cell = k_gauss_grad(dim,gamma)
%K_GAUSS_GRAD Summary of this function goes here
% 
% [kernel_grad_cell] = K_GAUSS_GRAD(gamma) 
% kernel_grad_cell = {grad_k_x, grad_k_y, ...} 
% each term is a func_handle
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: yihu $	$Date: 2016/01/21 21:03:30 $	$Revision: 0.1 $

kernel_grad_cell = cell(1,dim);
for i=1:dim
    kernel_grad_cell{i} = @(X1, X2) gradient_component(X1,X2,i,gamma);
end

end

function k_grad_x = gradient_component(X1, X2, i, gamma)

if length(X1)==size(X1,1)
  X1 = X1';
end;
if length(X2)==size(X2,1)
  X2 = X2';
end;

K = k_gauss(X1, X2, gamma);
differ = bsxfun(@minus, X1(i,:)', X2(i,:));
k_grad_x = -2*gamma*bsxfun(@times, K, differ);

end
