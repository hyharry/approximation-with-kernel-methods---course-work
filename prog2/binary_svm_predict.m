function [ytest, ftest] = binary_svm_predict(Ktest, y, alpha, b, tag)

if isrow(y)
    y = y(:);
end
ftest = Ktest*(alpha.*y) + b;

ytest = sign(ftest);

% for classification problem with values other than +1 and -1
% y = tag_1 if f > 0
% y = tag_2 if f other
if nargin>4
    n = length(ytest);
    for i=1:n
        if ytest(i)>0
            ytest(i) = tag(1);
        else
            ytest(i) = tag(2);
        end 
    end
end
