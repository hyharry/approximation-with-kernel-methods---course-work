function K = k_gauss(X1,X2,gamma)
%function K = k_gauss(X1,X2,gamma)
%
% function computing a Gaussian kernel matrix for the given data in
% the columns of X1, X2
%
% Bernard Haasdonk 17.10.2015

if length(X1)==prod(size(X1))
  X1 = X1(:);
end;
if length(X2)==prod(size(X2))
  X2 = X2(:);
end;

n1 = size(X1,2);
n2 = size(X2,2);

% we make use of the decomposition |x1-x2|^2 = |x1|^2 - 2 x1 x2 + |x2|^2
X1sqr = sum(X1.^2);
X2sqr = sum(X2.^2);
Dsqr = X1sqr(:) * ones(1,n2)  - 2 * X1' * X2 + ones(n1,1) *  X2sqr(:)';
K = exp(- gamma * Dsqr);