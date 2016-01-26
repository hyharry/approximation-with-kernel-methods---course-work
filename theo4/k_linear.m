function K = k_linear(X1,X2)
%function K = k_linear(X1,X2)
%  
% function computing a linear kernel matrix for the given data in
% the columns of X1, X2
%
% Bernard Haasdonk 18.11.2015

% if length(X1)==prod(size(X1))
%   X1 = X1(:);
% end;
% if length(X2)==prod(size(X2))
%   X2 = X2(:);
% end;

if length(X1)==size(X1,1)
  X1 = X1';
end;
if length(X2)==size(X2,1)
  X2 = X2';
end;

K = X1' * X2;