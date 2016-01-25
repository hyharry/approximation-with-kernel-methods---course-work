function [X_b, y_b] = multi_to_binary(X,y,i,j)
% get data set X_b according to class i and class j
%  transfomr y into binary form

I_c1 = find(y == i);
I_c2 = find(y == j);
X_c1 = X(:,I_c1);
X_c2 = X(:,I_c2);
X_b = [X_c1 X_c2];
y_b = [ones(1,length(I_c1)) -ones(1,length(I_c2))];