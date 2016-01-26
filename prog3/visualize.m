function visualize(f, Xvis, Yvis, Xcenter, gamma)
% f: approx function
% Xvis: visualization points
% Xcenter: center for f evaluation
% gamma: kernel parameter
% plot_prec: plot parameter

% $Author: yihu $	$Date: 2016/01/21 15:00:30 $	$Revision: 0.1 $

z = f([Xvis(:), Yvis(:)]);
Zvis = reshape(z,size(Xvis));

% Plot
surf(Xvis,Yvis,Zvis);
hold on
% scatter3(Xcenter(:,1), Xcenter(:,2),f(Xcenter),'r','filled');
scatter(Xcenter(:,1), Xcenter(:,2),'r','filled');

if gamma == 0
    title('exact');
else
    title(['\gamma =', num2str(gamma)]);
end

end