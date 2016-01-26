function [f_approx,C,beta,X_index, err_crit_sequence] = ...
    greedy_kernel_interpol(X,y,k,m_max,eps,mode)
%function [f_approx,C,beta,X_index, err_crit_sequence] = greedy_kernel_interpol(X,y,k,m_max,eps,mode)
%
% function doing greedy kernel interpolation
% mode is 'f', 'P' or 'f/P' for selection criterion

% B. Haasdonk 18.12.2015

y = y(:);
K = k(X,X);
Kxx = diag(K);

C = []; D = [];
beta = [];
sum_i_Ni_sqr_X = y*0;
X_index = [];
remainder = y;
not_selected = y*1;

m = 0;
err_crit = 1e10;
err_crit_sequence = [];
while (m < m_max) && (err_crit > eps)
  m = m+1;
  not_selected_index = find(not_selected);
  denominator = (Kxx - sum_i_Ni_sqr_X);  
  switch mode
   case 'f/P'
    % the following workaround switch to f-greedy for small P values:
    %    i = find(denominator<1e-8);
    %    denominator(i) = 1; 
    [err_crit,ind ] = max(remainder(not_selected_index).^2 .* ...
			  (denominator(not_selected_index).^(-1)));
   case 'f'
    [err_crit,ind ] = max(remainder(not_selected_index).^2);
   case 'P'
    [err_crit,ind ] = max(denominator(not_selected_index));
  end;
  x_m_index = not_selected_index(ind(1));
  err_crit_sequence = [err_crit_sequence, err_crit];
  
  %  if i==5
  %    keyboard;
  %  end;
  
  beta_m = (remainder(x_m_index))...
      /sqrt(Kxx(x_m_index)-sum_i_Ni_sqr_X(x_m_index));
  % extend coefficient vector for newton basis:
  beta = [beta; beta_m];
  
  % extend Newton basis 
  % a) update Cholesky factor D of K_m, i.e. D * D^T = K_m
  k_m = K(X_index, x_m_index);
  k_mm = Kxx(x_m_index);
  d_m = D \ k_m;
  d_mm = sqrt(k_mm - d_m' * d_m); 
  D = [D    , zeros(m-1,1); ...
       d_m' , d_mm         ];
  % b) update inverse Cholesky C Factor of K_m, i.e. C := D^-1:
  c_mm = d_mm^(-1);
  c_m = - c_mm*C'*d_m;
  C = [C    , zeros(m-1,1); ...
       c_m' , c_mm ];  
  
  % extend point index set:
  X_index = [X_index, x_m_index];
  not_selected(x_m_index) = 0;
  
  % evaluate new Newton basis vector in all points:
  Nm_X = K(:,X_index) * C(m,:)';
  sum_i_Ni_sqr_X = sum_i_Ni_sqr_X + Nm_X.^2;
  
  % update remainder
  remainder = remainder - Nm_X * beta_m;
  remainder(X_index) = 0; % set to 0 for not selection later
  
end;

f_approx = @(Xtest) f_hat(Xtest,X(:,X_index),C,beta,k);

function y_appr = f_hat(Xtest,Xm,C,beta,k)
Ktest = k(Xtest,Xm);
y_appr = zeros(size(Ktest,1),1);
for i = 1:length(beta)
  y_appr = y_appr + Ktest * C(i,:)'* beta(i);
end;
