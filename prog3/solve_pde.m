function f_approx = solve_pde(Xcenter, model, kernel, assemble_method)
% Solve pde, assemble_method for A,b assemble,
% construct approximant,
% 
% Xcenter: kernel center, one node per row
% model: poisson or generel
% assemble_method: RBF or other methods

% $Author: yihu $	$Date: 2016/01/21 21:56:44 $	$Revision: 0.1 $

%     kernel = @(X1,X2) k_gauss(X1,X2,gamma);
%     kernel_grad = k_gauss_grad(dim,gamma);
%     [A,b] = assemble_method(model,kernel,kernel_grad,Xcenter,Xcenter);
[A,b] = assemble_method(model,kernel,Xcenter,Xcenter);
alpha = A\b;
f_approx = @(Xtest) kernel(Xtest, Xcenter)*alpha;

end

