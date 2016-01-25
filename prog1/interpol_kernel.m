function f_approx = interpol_kernel( X, y, k, para )
%f_approx using Kernel Interpolation
%   input: X = [x1, x2, .. ,xn] ... n training points
%          y = [y1, y2, .. ,yn] ... n eval at n points
%          k in {l,p,d}         ... kernel method

if k=='l'
    para = 1;
end

K = kernel_matrix(X,X,k,para);
coeff = K\y';
f_approx = @(X_p) kernel_matrix(X_p,X,k,para)*coeff;

% switch k
%     case 'l'
%         K = kernel_matrix(X,X,'l');
%         coeff = K\y;
%         f_approx = @(X_p) kernel_matrix(X,X_p,'l');
%     case 'p'
%     case 'g'
%     otherwise
%         error('please specify a valid method (l,p,d)')
% end

end

