function res = interpolation_experiments(step)
%function res = interpolation_experiments(step)
%
% step 1 = spectral error decay for N_k and |inf| if f is in N_k, spaces nested
% step 2 = condition number degeneration with increasing n
% step 3 = reproduction of convergence order beta for kernel r^beta

% B. Haasdonk 16.12.2015

if nargin < 1
  step = 1;
end;

switch step
 case 1 % error decay for N_k and |inf| if f is in N_k, spaces nested
  gamma = 20;
  k = @(X1,X2)  k_gauss(X1,X2,gamma);
  xf = [0.123,0.965,0.5342];
  alphaf = [1;-2;4];
  f = @(x) k_gauss(x,xf,gamma)*alphaf;

  ns = [3,5,9,17,33,65,129];
  
  [err_linfty,err_Nk] = interpol_and_error(f,alphaf,xf,k,ns);

  plot(ns, [err_Nk;err_linfty]);
  legend('RKHS error norm','infty error norm');
  set(gca,'Yscale','log');
  set(gca,'Xscale','log');
  xlabel('number of train samples n');
  ylabel('error');
  
  disp('spectral convergence for gaussian kernel')
  disp('monotonicity for N_k error norm')
  disp('possible nonmonotonicity for l-inftyerror norm')
  disp('press enter to continue');
  pause

  figure;
  hs = 1./(ns-1);
  bound = 10 * sqrt(exp(- 10^-1 * abs (log(hs))./hs));
  plot(ns, [err_linfty;bound]);
  title('error and (qualitative) bound for Gaussian kernel')
  legend('infty error norm','sqrt(F(h))');
  set(gca,'Yscale','log');
  set(gca,'Xscale','log');
  xlabel('number of train samples n');
  ylabel('error');

  disp('linfty bound qualitatively scaling with bound')
 
 case 2 % condition number degeneration for settings of case 1
  
  gamma = 20;
  k = @(X1,X2)  k_gauss(X1,X2,gamma);
  ns = [3,5,9,17,33,65,129];
  cond_number = zeros(1,length(ns));
  for ni = 1:length(ns);
    n = ns(ni); 
    X = linspace(0,1,n); 
    K = k(X,X);
    cond_number(ni) = cond(K,2);
  end;  
  plot(ns,cond_number);
  title('condition number \kappa_2(K)');
  xlabel('number of train samples n');
  ylabel('\kappa_2(K)');
  set(gca,'Yscale','log');
  keyboard;
  disp('condition number terribly decaying')
  
 case 3 % reproduction of convergence order beta for kernel r^beta
  
  betas = [1,3,5];
%  betas = [2,3,4,5];
  err_linfty = cell(1,length(betas));
  for bi = 1:length(betas);
    beta = betas(bi);
    k = @(X1,X2)  k_neg_distance(X1,X2,beta);
    xf = [0.123,0.965,0.5342];
    alphaf = [1;-2;4];
    f = @(x) k(x,xf)*alphaf;
    
    ns = [3,5,9,17,33,65,129];
    
    err_linfty{bi}= interpol_and_error(f,alphaf,xf,k,ns);
  end;
  
  plot(ns, [err_linfty{1};err_linfty{2};err_linfty{3}]);
  title('linfty norm for r^\beta kernel')
  legend('\beta=1','\beta = 3','\beta = 5');
  set(gca,'Yscale','log');
  set(gca,'Xscale','log');
  xlabel('number of train samples n');
  ylabel('error');

  disp('convergence order for kernel r^beta')
  disp('order 1/2 for beta=1 satisfied')
  disp('order 3/2 for beta=3 oversatisfied')
  disp('order 5/2 for beta=5 oversatisfied')
  
 case 4 % problem of interpolation if f is not in RKHS?
  
  gamma = 20;
  k = @(X1,X2)  k_gauss(X1,X2,gamma);
%  xf = [0.123,0.965,0.5342];
%  alphaf = [1;-2;4];
%  f = @(x) k_gauss(x,xf,gamma)*alphaf;
  f = @(x) abs(x'-0.5);

  ns = [3,5,9,17,33,65,129];
  alphaf = [];
  xf = [];
  [err_linfty,err_Nk,fvals] = interpol_and_error(f,alphaf,xf,k,ns);

  plot(ns, [err_linfty]);
  legend('infty error norm');
  set(gca,'Yscale','log');
  set(gca,'Xscale','log');
  xlabel('number of train samples n');
  ylabel('error');
  
  figure;
  xtest = linspace(0,1,501);
  plot(xtest,fvals(1:6,:)');
  legend('f','n=3','n=5','n=9','n=17','n=33');
  xlabel('x coordinate');
  ylabel('function and interpolants');
  
  disp('problem if f_target is not in kernel space!');
  
 case 5 % greedy interpolation and convergence 
  
  gamma = 20;
  k = @(X1,X2)  k_gauss(X1,X2,gamma);
  xf = [0.123,0.965,0.5342];
  alphaf = [1;-2;4];
  f = @(x) k_gauss(x,xf,gamma)*alphaf;
  xf = [],
  alphaf = [];
  f = @(x) abs(x'-0.5);
%   f = @(x) sin(2*pi*x(:));
  
  n = 129;
  X = linspace(0,1,n); 
  y = f(X);
  Xtest = 0:1/500:1;
  f_test = f(Xtest);
%   greedy_mode = 'f/P';
 greedy_mode = 'f';
%  greedy_mode = 'P';
  figure;
  m_max = 20;
  m = 0;
  eps = 1e-10;
  err_crit_sequence = 1e10;
  while (m < m_max) & (err_crit_sequence(end)>eps) 
    m = m+1;
    [f_appr, C, beta,X_index,err_crit_sequence] = ...
	greedy_kernel_interpol(X,y,k,m,eps,greedy_mode);
    %  [XX_test,YY_test] = meshgrid(xx_test,xx_test);
    %  Xtest = [XX_test(:), YY_test(:)]';
    f_appr_test = f_appr(Xtest);
    %  fvals = [fvals;f_appr_test(:)'];
    %  err_linfty(ni) = max(abs(f_appr_test-f_test)); 
    err_linfty = max(abs(f_appr_test-f_test)); 
    cla;
    plot(Xtest',[f_appr_test,f_test]');
    hold on;
    plot(X(X_index),y(X_index),'ro');
    plot(X(X_index(end)),y(X_index(end)),'ro','Markerfacecolor',[1,0,0]);
    title(['iteration ',num2str(m)]);
    legend('f_m','f');
    disp(['maximum criterion before last extension: ',...
	  num2str(err_crit_sequence(end))]);
    disp('press enter to continue');
    pause;
  end;
  figure;
  plot(err_crit_sequence);
  set(gca,'Yscale','log');
  title('selection criterion development')
  xlabel('greedy iterations');
  %  err_Nk(ni) = sqrt(err_norm_sqr);
   
  % Show that f/P is optimal concerning N_k error convergence
  % compare with respect to linfty error

 case 6 % comparison of error decrease: greedy variants

  gamma = 20;
  k = @(X1,X2)  k_gauss(X1,X2,gamma);
  xf = [0.123,0.965,0.5342];
  alphaf = [1;-2;4];
  f = @(x) k_gauss(x,xf,gamma)*alphaf;
  %xf = [],
  %alphaf = [];
  %f = @(x) abs(x'-0.5);
  %f = @(x) sin(2*pi*x(:));
    
  n = 129;
  X = linspace(0,1,n); 
  y = f(X);
  Xtest = 0:1/500:1;
  f_test = f(Xtest);
  greedy_modes = {'f/P', 'f', 'P'};
  nmodes = length(greedy_modes);
  m_max = 20;
  eps = 1e-15;
  linfty_err_sequences = cell(1,length(nmodes));
  Nk_err_sequences = cell(1,length(nmodes));
  for mi = 1:nmodes;
    linfty_err_sequence = [];
    Nk_err_sequence = [];
    err_crit_sequence = 1e10;
    m = 0;
    greedy_mode = greedy_modes{mi};
    while (m < m_max) & (err_crit_sequence(end)>eps) 
      m = m+1;
      [f_appr, C, beta,X_index,err_crit_sequence] = ...
	  greedy_kernel_interpol(X,y,k,m,eps,greedy_mode);
      f_appr_test = f_appr(Xtest);
      linfty_err = max(abs(f_appr_test-f_test)); 
      linfty_err_sequence = [linfty_err_sequence,linfty_err];
      alpha = C' * beta;
      Nk_err = Nk_Omega_error(X(X_index),xf,k,alpha,alphaf);
      Nk_err_sequence = [Nk_err_sequence,Nk_err];      
    end;
    linfty_err_sequences{mi} = linfty_err_sequence;
    Nk_err_sequences{mi} = Nk_err_sequence;
  end;
  % plot results
  % linfty error
  figure;
  cols = 'bgrky';
  for i = 1:length(linfty_err_sequences)
    err_seq = linfty_err_sequences{i};
    plot(1:length(err_seq),err_seq,cols(i));
    hold on;
  end;
  set(gca,'Yscale','log');
  legend(greedy_modes);
  xlabel('greedy iteration');
  ylabel('error');
  title('l_\infty-error decay');

  % Nk error
  if ~isempty(alphaf)
    figure;
    cols = 'bgrky';
    for i = 1:length(Nk_err_sequences)
      err_seq = Nk_err_sequences{i};
      plot(1:length(err_seq),err_seq,cols(i));
      hold on;
    end;
    set(gca,'Yscale','log');
    legend(greedy_modes);
    xlabel('greedy iteration');
    ylabel('error');
    title('N_k(\Omega)-error decay');  
  end;

% case 7 % greedy convergence for q > 1
%
%  ds = [1,2];
%  gamma = 20;
%  k = @(X1,X2)  k_gauss(X1,X2,gamma);
%  fs = cell(1,length(ds));
%  alphaf = [1;-2;4];
%  linfty_err_sequences = cell(1,length(ds));
%  Nk_err_sequences = cell(1,length(ds));
%  
%  for di = 1:length(ds);
%    d = ds(di);
%    xf = ones(d,3)* diag([0.123,0.965,0.5342]);
%    f = @(x) k_gauss(x,xf,gamma)*alphaf;
%    %xf = [],
%    %alphaf = [];
%    %f = @(x) abs(x'-0.5);
%    %f = @(x) sin(2*pi*x(:));
%    
%    n = 5000;
%    %  X = linspace(0,1,n); 
%    X = rand(d,n);
%    y = f(X);
%    ntest = n;
%    Xtest = X;
%%    ntest = 500;
%%    Xtest = rand(d,ntest);
%    f_test = f(Xtest);
%    greedy_mode = 'f';
%    m_max = 20;
%    eps = 1e-15;
%    nmodes = 1;
%    linfty_err_sequence = [];
%    Nk_err_sequence = [];
%    err_crit_sequence = 1e10;
%    m = 0;
%    mi = 1;
%    while (m < m_max) & (err_crit_sequence(end)>eps) 
%      m = m+1;
%      [f_appr, C, beta,X_index,err_crit_sequence] = ...
%	  greedy_kernel_interpol(X,y,k,m,eps,greedy_mode);
%      f_appr_test = f_appr(Xtest);
%      linfty_err = max(abs(f_appr_test-f_test)); 
%      linfty_err_sequence = [linfty_err_sequence,linfty_err];
%      alpha = C' * beta;
%      Nk_err = Nk_Omega_error(X(:,X_index),xf,k,alpha,alphaf);
%      Nk_err_sequence = [Nk_err_sequence,Nk_err];      
%    end;
%    linfty_err_sequences{di} = linfty_err_sequence;
%    Nk_err_sequences{di} = Nk_err_sequence;
%  end;
%  % plot results
%  % linfty error
%  figure;
%  cols = 'bgrky';
%  for i = 1:length(linfty_err_sequences)
%    err_seq = linfty_err_sequences{i};
%    plot(1:length(err_seq),err_seq,cols(i));
%    hold on;
%  end;
%  set(gca,'Yscale','log');
%  legend(num2str(ds));
%  xlabel('greedy iteration');
%  ylabel('error');
%  title('l_\infty-error decay');
%
%  % Nk error
%  if ~isempty(alphaf)
%    figure;
%    cols = 'bgrky';
%    for i = 1:length(Nk_err_sequences)
%      err_seq = Nk_err_sequences{i};
%      plot(1:length(err_seq),err_seq,cols(i));
%      hold on;
%    end;
%    set(gca,'Yscale','log');
%    legend(num2str(ds));
%    xlabel('greedy iteration');
%    ylabel('error');
%    title('N_k(\Omega)-error decay');  
%  end;

  
 otherwise
  error('step unknown');
end;

function [err_linfty,err_Nk,fvals] = interpol_and_error(f,alphaf,xf,k,ns)
err_Nk = zeros(1,length(ns));
err_linfty = zeros(1,length(ns));
fvals = zeros(0,501);
for ni = 1:length(ns);
  n = ns(ni); 
  X = linspace(0,1,n); 
  y = f(X);
  [f_appr, alpha] = kernel_interpol(X,y,k);
  Xtest = 0:1/500:1;
  %  [XX_test,YY_test] = meshgrid(xx_test,xx_test);
  %  Xtest = [XX_test(:), YY_test(:)]';
  f_appr_test = f_appr(Xtest);
  fvals = [fvals;f_appr_test(:)'];
  f_test = f(Xtest);
  err_linfty(ni) = max(abs(f_appr_test-f_test)); 
  err_Nk(ni) = Nk_Omega_error(X,xf,k,alpha,alphaf);
end;
fvals = [f_test(:)';fvals];

function err = Nk_Omega_error(X,xf,k,alpha,alphaf)
x_err = [X,xf];
K = k(x_err,x_err);
alpha_err = [alpha;-alphaf];
err_norm_sqr = max(alpha_err' * K * alpha_err,0);
err = sqrt(err_norm_sqr);
