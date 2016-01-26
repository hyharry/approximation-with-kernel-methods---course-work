function error = post_process(f_exact, f_approx_list, gamma_list, Xvis, Yvis, Xcenter)
% Post Processing: error calculation and visualization
% 
% Xcenter: kernel center, one node per row
% Xvis: points for plotting
% gamma_list: different kernel parameters

% $Author: yihu $	$Date: 2016/01/21 21:56:44 $	$Revision: 0.1 $

test_number = length(gamma_list);
error = ones(1,test_number);

m = plot_arrange(test_number+1);
n = (test_number+1)/m;
for i=1:test_number
    f = f_approx_list{i};
    subplot(m,n,i);
    visualize(f,Xvis,Yvis,Xcenter,gamma_list(i));
    error(i) = max(abs(f_exact([Xvis(:),Yvis(:)]) - f([Xvis(:),Yvis(:)])));
end

% Plot exact solution
subplot(m,n,i+1);
visualize(f_exact,Xvis,Yvis,Xcenter,0);

% Plot error agains gamma
figure;
plot(gamma_list,error);
set(gca,'Xscale','log');
title('error of different gamma')
xlabel('gamma');
ylabel('L inf error')

% print errors
display(sprintf(repmat('-',1,70)));
for i = 1:test_number
    display(sprintf('gauss_kernel, gamma = %.3f, error = %5.4f', ...
        gamma_list(i), error(i)));
end

[best_gamma, ind] = min(error);
display(sprintf('best gamma = %.3f, error = %.4f', gamma_list(ind), best_gamma));

end

function m = plot_arrange(num)
mid = floor(sqrt(num));
while mod(num,mid) ~= 0
    mid = mid-1;
end
m = mid;
end