function kernel_laplace = k_gauss_laplace(gamma)
%K_GAUSS_GRAD 
% 
% [kernel_laplace] = K_GAUSS_GRAD(gamma)
% return is a function handel

% $Author: yihu $	$Date: 2016/01/21 21:03:30 $	$Revision: 0.1 $

term_1 = @(X1, X2) -2*gamma*k_gauss(X1,X2,gamma);

% Each row of X1 should be a point coord
Xdist = @(X1, X2) bsxfun(@plus,sum(X1.^2,2),sum(X2.^2,2)') - 2*X1*X2';

term_2 = @(X1, X2) 4*gamma^2*bsxfun(@times, k_gauss(X1,X2,gamma), Xdist(X1,X2));

kernel_laplace = @(X1, X2) term_1(X1, X2) + term_2(X1, X2);
