function [error] = solve_pde(Xcenter, Xvis, Yvis, model, kernel, assemble_method, gamma_list)
% Solve pde according to different gamma, use assemble_method for A,b assemble,
% display errors
% 
% Xcenter: kernel center, one node per row
% Xvis: points for plotting
% plot_prec: plot parameter
% model: poisson or generel
% gamma_list: different kernel parameters
% assemble_method: RBF or other methods
% label: label for assemble method

% $Author: yihu $	$Date: 2016/01/21 21:56:44 $	$Revision: 0.1 $

dim = size(Xcenter,2);

test_number = length(gamma_list);
error = ones(1,test_number);

for i=1:test_number
    gamma = gamma_list(i);
%     kernel = @(X1,X2) k_gauss(X1,X2,gamma);
%     kernel_grad = k_gauss_grad(dim,gamma);
%     [A,b] = assemble_method(model,kernel,kernel_grad,Xcenter,Xcenter);
    [A,b] = assemble_method(model,kernel,Xcenter,Xcenter);
    alpha = A\b;
    f = @(Xtest) kernel(Xtest, Xcenter)*alpha;
    visualize(f,Xvis,Yvis,Xcenter,gamma);
    error(i) = max(abs(model.solution([Xvis(:),Yvis(:)]) - f([Xvis(:),Yvis(:)])));
end

% Plot exact solution
% visualize(model.solution,Xvis,Yvis,Xcenter,0);

% print errors
display(sprintf(repmat('-',1,70)));
for i = 1:test_number
    display(sprintf('gauss_kernel, gamma = %.3f, error = %.4f', ...
        gamma_list(i), error(i)));
end

end

