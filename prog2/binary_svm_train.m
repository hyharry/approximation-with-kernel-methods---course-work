function [alpha, b] = binary_svm_train(K, y, C)

if isrow(y)
    y = y(:);
end
n = size(K,1);
H = (y*y').*K;
f = -ones(n,1);
Aeq = y';
beq = 0;
LB = zeros(n,1);
UB = C*ones(n,1);

alpha = quadprog(H,f,[],[],Aeq,beq,LB,UB);

b = y - K*(alpha.*y);

% take the averay of b vector
b = sum(b)/length(b);