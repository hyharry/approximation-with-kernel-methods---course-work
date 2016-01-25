function K = kernel_matrix( X, X_p, method, para)
% Generation of Kernel Matrix 
%   input: X, X_p .... two point sets
%                      [x1, x2, ... xn;
%                       y1, y2, ....yn;
%                       ..............]
%          method .... kernel
%          para   .... kernel parameters

if nargin == 2
    error('please enter the name of kernel method');
end

switch method
    case 'l'
        K = X'*X_p;
    case 'p'
        if nargin == 3
            error('please specify parameters (p,a)')
        end
        if length(para) == 1
            K = power(X'*X_p, para(1));
        else
            K = power(X'*X_p+para(2), para(1));
        end
    case 'g'
        if nargin == 3
            error('please specify parameters (gamma)')
        end
	% please mind that repmat is much slower than
	%  matrix multiplication
        X_abs = repmat(diag(X'*X), 1, size(X_p,2));
        X_p_abs = repmat(diag(X_p'*X_p)', size(X,2), 1);
        X_dot_X = X'*X_p;
        K = exp(-para * (X_abs+X_p_abs-2*X_dot_X));
    otherwise
            error('please enter a valid method (l,p,g)');
end
        
end

