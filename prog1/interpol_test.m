function error = interpol_test( m, k, para, m_p )
%INTERPOL_TEST: Test Interpolation for sin()sin() Function
%   input: m   ... equidistant points with distance 1/m
%          k   ... kernel method
%          para... kernel parameter
%          m_p ... test points to evaluate error
%   output:eror

% compute f for training points
n = m+1;
x = linspace(0.2,0.8,n);
[X1, X2] = meshgrid(x,x);
f = sin(2*pi*X1).*sin(4*pi*X2);

% transfer griddata to point coordinates
X = [X2(:)';X1(:)'];
y_tmp = permute(f,[2,1]); % y_tmp = f'
y = y_tmp(:);

% interpolation
% f_approx = interpol_kernel(X,y',k,para);
% row_mid_x = (X1(:,1:end-1) + X1(:,1:end_1))/2;
% row_mid_y = repmat(x',1,size(row_mid_x,2));

f_approx = interpol_kernel(X,transpose(y),k,para);

% plot interpolated result
xp = linspace(0,1,40);
[XP1,XP2] = meshgrid(xp,xp);
XP = [XP2(:)';XP1(:)'];
fp_approx_tmp = f_approx(XP);
fp_approx = reshape(fp_approx_tmp,40,40);
surf(XP1,XP2,fp_approx');

% fp = sin(2*pi*XP1).*sin(4*pi*XP2);
% 
% surf(XP1,XP2,fp)
% hold on
% scatter3(

error = 1;

end

