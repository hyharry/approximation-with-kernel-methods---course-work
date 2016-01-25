function Xtest = test_point_gen(a,b,n)
% generate equal dist grid for X_test

x = linspace(a,b,n);
[coord1, coord2] = meshgrid(x,x);
Xtest = [coord1(:)'; coord2(:)'];
